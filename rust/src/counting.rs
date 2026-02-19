use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Record};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

use crate::stats::fisher_strand_bias;
use crate::types::{BaseCounts, Variant};

use rayon::prelude::*;

use anyhow::{Context, Result};
use log::{debug, trace};
use bio::alignment::pairwise::Aligner;

/// Evidence accumulated for a single fragment (read pair) at a variant site.
/// Tracks the best base quality seen for each allele across both reads,
/// enabling quality-weighted consensus to resolve R1-vs-R2 conflicts.
///
/// Orientation tracking is per-allele: the strand direction stored for each
/// allele is the orientation of the read that provided the best-quality
/// evidence for that allele. This ensures Fisher's Strand Bias (FSB)
/// reflects the actual strand of the evidence, not just R1's strand.
#[derive(Debug, Clone)]
struct FragmentEvidence {
    /// Best base quality seen supporting REF across reads in this fragment
    best_ref_qual: u8,
    /// Best base quality seen supporting ALT across reads in this fragment
    best_alt_qual: u8,
    /// Orientation of the read providing best REF evidence
    best_ref_orientation: Option<bool>,
    /// Orientation of the read providing best ALT evidence
    best_alt_orientation: Option<bool>,
    /// Orientation of read 1 (fallback for orientation when allele-specific is unavailable)
    read1_orientation: Option<bool>,
    /// Orientation of read 2 (fallback)
    read2_orientation: Option<bool>,
}

impl FragmentEvidence {
    fn new() -> Self {
        FragmentEvidence {
            best_ref_qual: 0,
            best_alt_qual: 0,
            best_ref_orientation: None,
            best_alt_orientation: None,
            read1_orientation: None,
            read2_orientation: None,
        }
    }

    /// Record a read's allele call and orientation into this fragment's evidence.
    ///
    /// Tracks per-allele orientation: when a new best-quality observation is
    /// recorded for REF or ALT, the orientation of THAT read is stored.
    /// This couples the strand direction to the winning evidence, not just R1.
    fn observe(&mut self, is_ref: bool, is_alt: bool, base_qual: u8, is_read1: bool, is_forward: bool) {
        if is_ref && base_qual > self.best_ref_qual {
            self.best_ref_qual = base_qual;
            self.best_ref_orientation = Some(is_forward);
        }
        if is_alt && base_qual > self.best_alt_qual {
            self.best_alt_qual = base_qual;
            self.best_alt_orientation = Some(is_forward);
        }
        // Track R1/R2 orientation as fallback
        if is_read1 {
            self.read1_orientation = Some(is_forward);
        } else {
            self.read2_orientation = Some(is_forward);
        }
    }

    /// Resolve this fragment's allele call using quality-weighted consensus.
    /// Returns (is_ref, is_alt) — exactly one will be true.
    fn resolve(&self, qual_diff_threshold: u8) -> (bool, bool) {
        let has_ref = self.best_ref_qual > 0;
        let has_alt = self.best_alt_qual > 0;

        match (has_ref, has_alt) {
            (true, false) => (true, false),   // Only REF evidence
            (false, true) => (false, true),   // Only ALT evidence
            (true, true) => {
                // Conflict: both alleles seen across reads in this fragment.
                // Use quality-weighted consensus: higher quality wins.
                // If quality difference is within threshold, discard the fragment
                // to avoid biasing VAF in either direction — critical for low-VAF
                // cfDNA detection where every fragment matters.
                if self.best_ref_qual > self.best_alt_qual + qual_diff_threshold {
                    (true, false)  // REF wins by quality margin
                } else if self.best_alt_qual > self.best_ref_qual + qual_diff_threshold {
                    (false, true)  // ALT wins by quality margin
                } else {
                    // Within threshold — ambiguous, discard to preserve VAF accuracy
                    (false, false)
                }
            }
            (false, false) => (false, false),  // Should not happen (filtered earlier)
        }
    }

    /// Get orientation for REF allele. Uses the strand of the read that provided
    /// the best REF evidence, falling back to R1 > R2 if not available.
    fn ref_orientation(&self) -> Option<bool> {
        self.best_ref_orientation
            .or(self.read1_orientation)
            .or(self.read2_orientation)
    }

    /// Get orientation for ALT allele. Uses the strand of the read that provided
    /// the best ALT evidence, falling back to R1 > R2 if not available.
    fn alt_orientation(&self) -> Option<bool> {
        self.best_alt_orientation
            .or(self.read1_orientation)
            .or(self.read2_orientation)
    }
}

/// Hash a QNAME to u64 for memory-efficient fragment tracking.
/// Using DefaultHasher for speed — collision probability is negligible
/// for typical variant-level read counts (~1000 fragments).
#[inline]
fn hash_qname(qname: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Count bases for a list of variants in a BAM file.
///
/// When `decomposed` is provided (same length as `variants`), variants with
/// a `Some(decomposed_variant)` are counted twice — once with the original
/// allele and once with the corrected allele. The result with the higher
/// `ad` (alt_count) is returned, with `used_decomposed` set accordingly.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10))]
pub fn count_bam(
    py: Python<'_>,
    bam_path: String,
    variants: Vec<Variant>,
    decomposed: Vec<Option<Variant>>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: usize,
    fragment_qual_threshold: u8,
) -> PyResult<Vec<BaseCounts>> {
    // We cannot share a single IndexedReader across threads because it's not Sync.
    // Instead, we use rayon's map_init to initialize a reader for each thread.
    // This is efficient because map_init reuses the thread-local state (the reader)
    // for multiple items processed by that thread.

    // Configure thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to build thread pool: {}", e)))?;

    // Zip variants with their decomposed counterparts for parallel iteration
    let paired: Vec<_> = variants.into_iter().zip(decomposed.into_iter()).collect();

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<BaseCounts>, anyhow::Error> = py.allow_threads(move || {
        pool.install(|| {
            paired
                .par_iter()
                .map_init(
                    || {
                        // Initialize thread-local BAM reader
                        bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                            anyhow::anyhow!("Failed to open BAM: {}", e)
                        })
                    },
                    |bam_result, (variant, decomp_opt)| {
                        // Get the reader or return error if initialization failed
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        let counts_orig = count_single_variant(
                            bam,
                            variant,
                            min_mapq,
                            min_baseq,
                            filter_duplicates,
                            filter_secondary,
                            filter_supplementary,
                            filter_qc_failed,
                            filter_improper_pair,
                            filter_indel,
                            fragment_qual_threshold,
                        )?;

                        // Dual-count: if a decomposed variant exists, count it too
                        // and return whichever has the higher alt_count.
                        if let Some(decomp) = decomp_opt {
                            let counts_decomp = count_single_variant(
                                bam,
                                decomp,
                                min_mapq,
                                min_baseq,
                                filter_duplicates,
                                filter_secondary,
                                filter_supplementary,
                                filter_qc_failed,
                                filter_improper_pair,
                                filter_indel,
                                fragment_qual_threshold,
                            )?;

                            if counts_decomp.ad > counts_orig.ad {
                                // Sanity: both hypotheses count the same reads at the
                                // same locus, so DP should be nearly identical. A large
                                // divergence indicates a counting bug.
                                debug_assert!(
                                    (counts_decomp.dp as i64 - counts_orig.dp as i64).abs() <= 2,
                                    "DP mismatch in dual-counting: decomp={} orig={} at {}:{} {}→{}",
                                    counts_decomp.dp, counts_orig.dp,
                                    variant.chrom, variant.pos + 1,
                                    variant.ref_allele, decomp.alt_allele
                                );
                                debug!(
                                    "Homopolymer decomp: corrected allele wins \
                                     (ad={} vs orig ad={}, dp_decomp={}, dp_orig={}) \
                                     for {}:{} {}→{}",
                                    counts_decomp.ad, counts_orig.ad,
                                    counts_decomp.dp, counts_orig.dp,
                                    variant.chrom, variant.pos + 1,
                                    variant.ref_allele, decomp.alt_allele,
                                );
                                return Ok(BaseCounts {
                                    used_decomposed: true,
                                    ..counts_decomp
                                });
                            }
                        }

                        Ok(counts_orig)
                    },
                )
                .collect()
        })
    });

    // Map anyhow::Error back to PyErr
    match results {
        Ok(r) => Ok(r),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e))),
    }
}

#[allow(clippy::too_many_arguments)]
fn count_single_variant(
    bam: &mut bam::IndexedReader,
    variant: &Variant,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    fragment_qual_threshold: u8,
) -> Result<BaseCounts> {
    let tid = bam.header().tid(variant.chrom.as_bytes()).ok_or_else(|| {
        anyhow::anyhow!("Chromosome not found in BAM: {}", variant.chrom)
    })?;

    // Fetch region around the variant. For windowed indel detection (±5bp),
    // we expand the window so that reads with shifted indels are also retrieved.
    let window_pad: i64 = 5;
    let start = (variant.pos - window_pad).max(0);
    let end = variant.pos + (variant.ref_allele.len() as i64) + window_pad;

    bam.fetch((tid, start, end)).context("Failed to fetch region")?;

    let mut counts = BaseCounts::default();

    // Fragment tracking: QNAME hash -> FragmentEvidence
    // Using u64 hash keys instead of String for memory efficiency.
    let mut fragments: HashMap<u64, FragmentEvidence> = HashMap::new();

    // Quality threshold for fragment consensus tiebreaking.
    // When R1 and R2 disagree, the allele with higher quality wins
    // only if the quality difference exceeds this threshold.
    // Configurable via --fragment-qual-threshold (default: 10).
    let qual_diff_threshold: u8 = fragment_qual_threshold;

    // Create SW aligners ONCE per variant, not per read (indelpost pattern).
    // bio::alignment::pairwise::Aligner reuses internal DP buffers on
    // subsequent calls, avoiding repeated O(n×m) heap allocation.
    let score_fn = |a: u8, b: u8| -> i32 {
        if a == b'N' || b == b'N' { 0 } else if a == b { 1 } else { -1 }
    };
    // ALT + REF: Same affine gap penalties for fair comparison.
    // With raw read window extraction (extract_raw_read_window), the read
    // includes insertion bases making it potentially longer than either
    // haplotype. Both aligners need gap tolerance for correct scoring.
    let mut alt_aligner = Aligner::new(-5, -1, &score_fn);
    let mut ref_aligner = Aligner::new(-5, -1, &score_fn);

    for result in bam.records() {
        let record = result.context("Error reading BAM record")?;

        // Apply filters
        if filter_duplicates && record.is_duplicate() {
            continue;
        }
        if filter_secondary && record.is_secondary() {
            continue;
        }
        if filter_supplementary && record.is_supplementary() {
            continue;
        }
        if filter_qc_failed && record.is_quality_check_failed() {
            continue;
        }
        if filter_improper_pair && !record.is_proper_pair() {
            continue;
        }

        // Indel filter: check CIGAR for Ins or Del
        if filter_indel {
            let has_indel = record.cigar().iter().any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)));
            if has_indel {
                continue;
            }
        }

        if record.mapq() < min_mapq {
            continue;
        }

        // Determine allele status and base quality at variant position
        let (is_ref, is_alt, base_qual) = check_allele_with_qual(
            &record, variant, min_baseq, &mut alt_aligner, &mut ref_aligner,
        );

        if !is_ref && !is_alt {
            continue;
        }

        // Update per-read counts (unchanged — these are read-level, not fragment-level)
        counts.dp += 1;
        let is_reverse = record.is_reverse();

        if is_reverse {
            counts.dp_rev += 1;
        } else {
            counts.dp_fwd += 1;
        }

        if is_ref {
            counts.rd += 1;
            if is_reverse {
                counts.rd_rev += 1;
            } else {
                counts.rd_fwd += 1;
            }
        } else if is_alt {
            counts.ad += 1;
            if is_reverse {
                counts.ad_rev += 1;
            } else {
                counts.ad_fwd += 1;
            }
        }

        // Accumulate fragment evidence using QNAME hash
        let qname_hash = hash_qname(record.qname());
        let is_read1 = record.is_first_in_template();
        let is_forward = !is_reverse;

        let evidence = fragments.entry(qname_hash).or_insert_with(FragmentEvidence::new);
        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward);
    }

    // Resolve fragment-level counts using quality-weighted consensus.
    // Each fragment contributes exactly ONE allele call (REF xor ALT),
    // preventing the double-counting bug where R1=REF + R2=ALT
    // inflated both rdf and adf.
    //
    // Strand bias uses allele-specific orientation: the strand of the read
    // that provided the best evidence for the winning allele, not just R1.
    // Example: R1=Fwd/REF(Q10) + R2=Rev/ALT(Q30) → ALT wins, counted as
    // adf_rev (not adf_fwd).
    for evidence in fragments.values() {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);

        // Count every fragment in dpf regardless of consensus outcome.
        // Discarded fragments (ambiguous R1-vs-R2 within quality threshold)
        // are still real molecules — tracking them in dpf makes the gap
        // dpf - (rdf + adf) a useful quality metric for the locus.
        counts.dpf += 1;

        if frag_ref {
            counts.rdf += 1;
            // Use REF-specific orientation (strand of best REF evidence)
            if let Some(ori) = evidence.ref_orientation() {
                if ori {
                    counts.rdf_fwd += 1;
                } else {
                    counts.rdf_rev += 1;
                }
            }
        } else if frag_alt {
            counts.adf += 1;
            // Use ALT-specific orientation (strand of best ALT evidence)
            if let Some(ori) = evidence.alt_orientation() {
                if ori {
                    counts.adf_fwd += 1;
                } else {
                    counts.adf_rev += 1;
                }
            }
        }
    }

    // Calculate stats
    let (sb_pval, sb_or) =
        fisher_strand_bias(counts.rd_fwd, counts.rd_rev, counts.ad_fwd, counts.ad_rev);
    counts.sb_pval = sb_pval;
    counts.sb_or = sb_or;

    let (fsb_pval, fsb_or) = fisher_strand_bias(
        counts.rdf_fwd,
        counts.rdf_rev,
        counts.adf_fwd,
        counts.adf_rev,
    );
    counts.fsb_pval = fsb_pval;
    counts.fsb_or = fsb_or;

    Ok(counts)
}

/// Check if a read supports the reference or alternate allele.
/// Returns (is_ref, is_alt, base_quality) where base_quality is the
/// quality score at the variant position (used for fragment consensus).
///
/// Each variant-type handler returns quality directly from its own CIGAR
/// walk, ensuring correct quality extraction even for reads carrying
/// indels at the variant position.
///
/// The `alt_aligner` and `ref_aligner` are reusable SW aligners created
/// once per variant in `count_single_variant()` and threaded through to
/// avoid per-read allocation (indelpost pattern).
fn check_allele_with_qual<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> (bool, bool, u8) {
    let variant_type = &variant.variant_type;
    debug!("check_allele type={} pos={} ref={} alt={}", variant_type, variant.pos, variant.ref_allele, variant.alt_allele);

    match variant_type.as_str() {
        "SNP" => check_snp(record, variant, min_baseq),
        "INSERTION" => check_insertion(record, variant, min_baseq, alt_aligner, ref_aligner),
        "DELETION" => check_deletion(record, variant, min_baseq, alt_aligner, ref_aligner),
        "MNP" => {
            // MNP strict match first; fall back to quality-masked haplotype
            // comparison (check_complex → Phase 3 SW) if inconclusive.
            // Without this fallback, a single low-quality base in the MNP
            // block would hard-reject the read, reducing sensitivity in
            // noisy samples (ctDNA, FFPE).
            let (is_ref, is_alt, qual) = check_mnp(record, variant, min_baseq);
            if !is_ref && !is_alt {
                debug!("check_mnp inconclusive, falling back to check_complex");
                check_complex(record, variant, min_baseq, alt_aligner, ref_aligner)
            } else {
                (is_ref, is_alt, qual)
            }
        },
        "COMPLEX" => check_complex(record, variant, min_baseq, alt_aligner, ref_aligner),
        _ => {
            // Auto-detect MNP if type is not explicit but looks like one
            if variant.ref_allele.len() == variant.alt_allele.len()
                && variant.ref_allele.len() > 1
            {
                debug!("Auto-detected MNP");
                let (is_ref, is_alt, qual) = check_mnp(record, variant, min_baseq);
                if !is_ref && !is_alt {
                    debug!("Auto-detected MNP inconclusive, falling back to check_complex");
                    check_complex(record, variant, min_baseq, alt_aligner, ref_aligner)
                } else {
                    (is_ref, is_alt, qual)
                }
            } else {
                // Fallback to complex check for anything else (e.g. DelIns)
                check_complex(record, variant, min_baseq, alt_aligner, ref_aligner)
            }
        }
    }
}




/// Extract contiguous **raw** read bases spanning a genomic window `[win_start, win_end)`.
///
/// Unlike `extract_read_subsequence` which concatenates CIGAR-projected bases
/// (breaking complex variants represented as DEL+INS), this function finds the
/// read position range overlapping the window and returns the **original read
/// sequence** — preserving the true underlying allele for SW alignment.
///
/// ## Why this matters for complex variants
///
/// For EPHA7 (REF=TCC, ALT=CT), BWA represents ALT reads as `91M2D6M4I`.
/// CIGAR-projected extraction produces `GGAAACTCTCCAAAA` which matches
/// neither REF (`GGAAATCCACTCC`) nor ALT (`GGAAACTACTCC`).
/// Raw extraction produces `GGAAACTCTCCAAAA` from the same read positions,
/// which SW alignment can correctly classify by scoring against both haplotypes.
///
/// ## Algorithm
///
/// 1. Walk CIGAR to find the first read_pos where ref enters the window
/// 2. Walk CIGAR to find the last read_pos where ref is still in the window
///    (including insertion bases at window boundaries)
/// 3. Return `seq[first_read_pos..=last_read_pos]`
///
/// Returns `(bases, quals)` — empty if the read doesn't overlap the window.
///
/// `variant_pos` and `variant_ref_len` restrict soft-clip inclusion to clips
/// adjacent to the variant site, preventing adapter/chimeric garbage from
/// large leading clips from polluting Phase 3 alignment.
fn extract_raw_read_window(
    record: &Record,
    win_start: i64,
    win_end: i64,
    variant_pos: i64,
    variant_ref_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;
    let seq = record.seq();
    let quals = record.qual();

    let mut first_read_pos: Option<usize> = None;
    let mut last_read_pos: Option<usize> = None;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let len_usize = *len as usize;

                // Check overlap with [win_start, win_end)
                let overlap_start = std::cmp::max(ref_pos, win_start);
                let overlap_end = std::cmp::min(ref_pos + len_i64, win_end);
                if overlap_start < overlap_end {
                    let offset = (overlap_start - ref_pos) as usize;
                    let count = (overlap_end - overlap_start) as usize;
                    let rp_start = read_pos + offset;
                    let rp_end = rp_start + count - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(rp_start);
                    }
                    last_read_pos = Some(rp_end);
                }
                ref_pos += len_i64;
                read_pos += len_usize;
            }
            Cigar::Ins(len) => {
                let len_usize = *len as usize;
                // Include insertion bases if they fall within the window boundary.
                // An insertion at ref_pos X appears between ref X-1 and X.
                if ref_pos >= win_start && ref_pos <= win_end {
                    let rp_end = read_pos + len_usize - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(read_pos);
                    }
                    last_read_pos = Some(rp_end);
                }
                read_pos += len_usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                let len_usize = *len as usize;
                // Only include soft clips adjacent to the variant position
                // itself (within ±1bp), not the broad fetching window.
                // Large leading clips (e.g. 100S) at the window boundary
                // would inject garbage adapter/chimeric sequence into Phase 3
                // SW alignment, dragging down scores and causing valid reads
                // to fail the alt_score >= ref_score + margin check.
                let var_end = variant_pos + variant_ref_len as i64;
                if ref_pos >= variant_pos - 1 && ref_pos <= var_end + 1 {
                    let rp_end = read_pos + len_usize - 1;
                    if first_read_pos.is_none() {
                        first_read_pos = Some(read_pos);
                    }
                    last_read_pos = Some(rp_end);
                }
                read_pos += len_usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // Extract contiguous raw read bases from first to last overlapping position
    match (first_read_pos, last_read_pos) {
        (Some(first), Some(last)) if first <= last && last < seq.len() => {
            let bases: Vec<u8> = (first..=last).map(|i| seq[i]).collect();
            let base_quals: Vec<u8> = quals[first..=last].to_vec();
            (bases, base_quals)
        }
        _ => (Vec::new(), Vec::new()),
    }
}

/// Quick pre-filter before Smith-Waterman alignment (indelpost pattern).
///
/// Returns `true` only if the read shows evidence that could change the
/// allele call via SW realignment:
/// - Soft-clipping within or near the variant window
/// - CIGAR contains indels within `[win_start, win_end)`
///
/// Clean M-only reads that span the variant are skipped (the vast majority
/// at any locus), eliminating ~80-90% of unnecessary SW calls.
fn is_worth_realignment(record: &Record, win_start: i64, win_end: i64) -> bool {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Ins(_) => {
                // Insertion at this ref_pos — check if within window
                if ref_pos >= win_start && ref_pos <= win_end {
                    return true;
                }
            }
            Cigar::Del(len) => {
                let del_end = ref_pos + *len as i64;
                // Deletion overlapping window
                if ref_pos < win_end && del_end > win_start {
                    return true;
                }
                ref_pos = del_end;
            }
            Cigar::SoftClip(_) => {
                // Soft-clip near the variant window suggests misalignment
                if ref_pos >= win_start - 5 && ref_pos <= win_end + 5 {
                    return true;
                }
            }
            Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            _ => {}
        }
    }
    false
}

/// Classify a read subsequence as REF or ALT using dual-haplotype
/// Smith-Waterman alignment (inspired by indelpost).
///
/// Builds ref/alt haplotypes from `variant.ref_context`, then aligns
/// `read_seq` against both using the provided reusable `Aligner` instances.
///
/// **Alignment mode**: Semiglobal — haplotype (query) must be fully aligned,
/// read (text) has free overhangs. This answers "does this read *contain*
/// this haplotype?" — the standard approach for pattern-in-text problems.
///
/// **Haplotype trimming**: When a haplotype exceeds the read length (common
/// with large `context_padding` and complex CIGARs), the haplotype is
/// symmetrically trimmed from flanks to fit. This ensures semiglobal alignment
/// penalizes only allele differences, not excess length. A fixed score margin
/// of 2 works reliably because shared flanking bases contribute equally to
/// both allele scores — trimming preserves this invariant.
///
/// **Memory optimization**: Aligners are created once per variant in
/// `count_single_variant()` and reused for all reads, avoiding repeated
/// O(n×m) DP matrix allocation (indelpost pattern).
///
/// Returns `(is_ref, is_alt, base_qual)` where `base_qual` is the median
/// quality across the read subsequence bases (GATK standard).
///
/// # Alignment semantics
/// Uses semiglobal(read, haplotype): the read is the pattern (globally consumed,
/// no free overhangs) and the haplotype is the text (free leading/trailing gaps).
/// This allows the aligner to find the best-scoring window within a longer
/// haplotype that matches the read — critical for large deletions where the
/// REF haplotype (e.g. 501bp for a 500bp deletion) far exceeds the read length
/// (~150bp), and for off-center reads that don't symmetrically overlap the variant.
fn classify_by_alignment<F: Fn(u8, u8) -> i32>(
    read_seq: &[u8],
    read_quals: &[u8],
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> (bool, bool, u8) {
    let ref_context = match &variant.ref_context {
        Some(ctx) => ctx.as_bytes(),
        None => return (false, false, 0),
    };

    // Build haplotypes from ref_context.
    // ref_context covers [ref_context_start, ref_context_start + len).
    // The variant sits at variant.pos within this context.
    let offset = (variant.pos - variant.ref_context_start) as usize;
    let ref_len = variant.ref_allele.len();

    // Guard: ensure offset is valid within ref_context
    if offset + ref_len > ref_context.len() {
        debug!(
            "classify_by_alignment: offset {} + ref_len {} exceeds context len {}",
            offset, ref_len, ref_context.len()
        );
        return (false, false, 0);
    }

    let left_ctx = &ref_context[..offset];
    let right_ctx = &ref_context[offset + ref_len..];

    // ref_hap = left_ctx + REF + right_ctx (should equal ref_context)
    let ref_hap: Vec<u8> = left_ctx
        .iter()
        .chain(variant.ref_allele.as_bytes())
        .chain(right_ctx.iter())
        .copied()
        .collect();

    // alt_hap = left_ctx + ALT + right_ctx
    let alt_hap: Vec<u8> = left_ctx
        .iter()
        .chain(variant.alt_allele.as_bytes())
        .chain(right_ctx.iter())
        .copied()
        .collect();

    // Mask low-quality bases as N so they don't bias scoring.
    let masked_seq: Vec<u8> = read_seq
        .iter()
        .zip(read_quals.iter())
        .map(|(&b, &q)| if q >= min_baseq { b } else { b'N' })
        .collect();

    // Skip if too few usable bases.
    let read_len = read_seq.len();
    let usable_count = read_quals.iter().filter(|&&q| q >= min_baseq).count();
    if usable_count < 3 {
        debug!("classify_by_alignment: only {} usable bases — skipping", usable_count);
        return (false, false, 0);
    }

    // Use the provided reusable aligners (created once per variant).
    trace!(
        "classify_by_alignment: read_len={} ref_hap_len={} alt_hap_len={} usable={}",
        read_len, ref_hap.len(), alt_hap.len(), usable_count
    );
    trace!(
        "classify_by_alignment seqs: read={} alt_hap={} ref_hap={}",
        String::from_utf8_lossy(&masked_seq),
        String::from_utf8_lossy(&alt_hap),
        String::from_utf8_lossy(&ref_hap)
    );

    // Semiglobal alignment: read = pattern (fully consumed, every base must
    // participate), haplotype = text (free leading/trailing overhangs).
    // The aligner slides the read along the haplotype to find the best window.
    // This correctly handles:
    //  - Large deletions: 150bp read vs 501bp REF haplotype
    //  - Off-center reads: no symmetric trimming needed
    //  - Short ALT haplotypes: read may overhang, but all read bases align
    let alt_aln = alt_aligner.semiglobal(&masked_seq, &alt_hap);
    let ref_aln = ref_aligner.semiglobal(&masked_seq, &ref_hap);

    let med_qual = median_qual(read_quals, min_baseq);

    // Score margin: require at least 2 points difference to call.
    // Uses >= (not >) because with larger context_padding, the SW aligner
    // can find shifted alignments that reduce the score difference from 3
    // to exactly 2. Strict > would incorrectly classify these as ambiguous.
    let margin = 2;
    let is_alt = alt_aln.score >= ref_aln.score + margin;
    let is_ref = ref_aln.score >= alt_aln.score + margin;

    debug!(
        "classify_by_alignment: alt_score={} ref_score={} margin={} is_ref={} is_alt={} \
         read_len={} ref_hap={} alt_hap={}",
        alt_aln.score, ref_aln.score, margin, is_ref, is_alt,
        read_len, ref_hap.len(), alt_hap.len()
    );

    if is_alt {
        (false, true, med_qual)
    } else if is_ref {
        (true, false, med_qual)
    } else {
        // Ambiguous — scores too close
        (false, false, 0)
    }
}

/// Check if a read supports a complex variant (indel + substitution).
///
/// Uses **haplotype reconstruction**: walks the CIGAR to rebuild what the read
/// shows for the genomic region covered by REF, then compares the reconstructed
/// sequence to both REF and ALT using **quality-aware masked comparison**.
///
/// ## Masked Comparison ("Reliable Intersection")
///
/// Instead of requiring exact byte-for-byte match, bases with quality below
/// `min_baseq` are **masked out** — they cannot vote for either allele. Only
/// "reliable" (high-quality) bases participate in the comparison.
///
/// Three cases based on reconstructed sequence length:
/// - **Case A** (`recon == alt == ref` length): simultaneous REF/ALT check with
///   ambiguity detection. If reliable bases match *both*, read is discarded.
/// - **Case B** (`recon == alt` length only): masked comparison against ALT only.
/// - **Case C** (`recon == ref` length only): masked comparison against REF only.
///
/// Returns (is_ref, is_alt, base_qual) where base_qual is the median quality
/// across the reconstructed haplotype bases, used for fragment consensus.
fn check_complex<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> (bool, bool, u8) {
    let start_pos = variant.pos;
    let end_pos = variant.pos + variant.ref_allele.len() as i64; // exclusive

    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;
    let seq = record.seq();
    let quals = record.qual();

    // Pre-allocate reconstruction buffers for performance
    let capacity = seq.len();
    let mut reconstructed_seq: Vec<u8> = Vec::with_capacity(capacity);
    let mut quals_per_base: Vec<u8> = Vec::with_capacity(capacity);

    debug!(
        "check_complex start: pos={} ref={} alt={}",
        start_pos, variant.ref_allele, variant.alt_allele
    );

    // --- Phase 1: Haplotype Reconstruction ---
    // Walk the CIGAR to reconstruct what the read shows for [start_pos, end_pos).
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let len_usize = *len as usize;

                // Intersection of [ref_pos, ref_pos + len) and [start_pos, end_pos)
                let overlap_start = std::cmp::max(ref_pos, start_pos);
                let overlap_end = std::cmp::min(ref_pos + len_i64, end_pos);

                if overlap_start < overlap_end {
                    let offset_in_op = (overlap_start - ref_pos) as usize;
                    let overlap_len = (overlap_end - overlap_start) as usize;
                    let current_read_pos = read_pos + offset_in_op;

                    for i in 0..overlap_len {
                        let p = current_read_pos + i;
                        if p >= seq.len() {
                            break;
                        }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                ref_pos += len_i64;
                read_pos += len_usize;
            }
            Cigar::Ins(len) => {
                let len_usize = *len as usize;
                if ref_pos >= start_pos && ref_pos <= end_pos {
                    for i in 0..len_usize {
                        let p = read_pos + i;
                        if p >= seq.len() {
                            break;
                        }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                read_pos += len_usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                let len_usize = *len as usize;
                // P1-2: Include soft-clipped bases that overlap the variant window.
                // Soft clips don't consume reference, so ref_pos is unchanged.
                // This recovers evidence from reads where the aligner clipped
                // the variant-supporting bases (inspired by VarDict's approach).
                if ref_pos >= start_pos && ref_pos < end_pos {
                    for i in 0..len_usize {
                        let p = read_pos + i;
                        if p >= seq.len() { break; }
                        reconstructed_seq.push(seq[p]);
                        quals_per_base.push(quals[p]);
                    }
                }
                read_pos += len_usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let reconstructed_str = String::from_utf8_lossy(&reconstructed_seq);
    debug!(
        "Reconstructed: '{}' (len={})",
        reconstructed_str,
        reconstructed_seq.len()
    );

    // Median quality across the reconstructed haplotype for fragment consensus.
    // This is the quality we return for whichever allele matches.
    let med_haplotype_qual = median_qual(&quals_per_base, min_baseq);

    // --- Phase 2: Quality-Aware Masked Comparison ---
    // Mask out low-quality bases. Only reliable bases (qual >= min_baseq) vote.
    let alt_bytes = variant.alt_allele.as_bytes();
    let ref_bytes = variant.ref_allele.as_bytes();
    let recon_len = reconstructed_seq.len();
    let matches_alt_len = recon_len == alt_bytes.len();
    let matches_ref_len = recon_len == ref_bytes.len();

    // Guard: pathologically short reconstruction for large-REF variants.
    // When a large region (e.g. 1024bp deletion) produces a tiny reconstruction
    // (e.g. 1bp), direct Phase 2 comparison is unreliable — a 1bp recon would
    // trivially match a 1bp ALT allele, causing overcounting.
    // Skip to Phase 3 (Smith-Waterman alignment) which uses full haplotype context.
    let ref_len = ref_bytes.len();
    if ref_len > 50 && recon_len > 0 && recon_len < ref_len / 10 {
        debug!(
            "check_complex: recon_len={} is <10% of ref_len={} — \
             skipping Phase 2 (unreliable direct comparison)",
            recon_len, ref_len
        );
    } else if matches_alt_len && matches_ref_len {
        // Case A: Equal-length REF and ALT — need simultaneous check + ambiguity detection
        let (mismatches_alt, mismatches_ref, reliable_count) =
            masked_dual_compare(&reconstructed_seq, &quals_per_base, alt_bytes, ref_bytes, min_baseq);

        debug!(
            "Case A: reliable={} mm_alt={} mm_ref={}",
            reliable_count, mismatches_alt, mismatches_ref
        );

        // Step 1: No reliable data → discard (MUST come first)
        if reliable_count == 0 {
            debug!("No reliable bases — discarding");
            return (false, false, 0);
        }

        // Step 2: Ambiguity — reliable bases match both alleles → discard
        if mismatches_alt == 0 && mismatches_ref == 0 {
            debug!("Ambiguous: reliable bases match both REF and ALT — discarding");
            return (false, false, 0);
        }

        // Step 3: Unambiguous match
        if mismatches_alt == 0 {
            debug!("Matches ALT on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return (false, true, med_haplotype_qual);
        }
        if mismatches_ref == 0 {
            debug!("Matches REF on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return (true, false, med_haplotype_qual);
        }

        // Step 4: Neither matches on reliable bases
        debug!("No match: mm_alt={} mm_ref={}", mismatches_alt, mismatches_ref);
    } else if matches_alt_len {
        // Case B: Only ALT length matches (e.g., DelIns) — no ambiguity possible
        let (mismatches, reliable_count) =
            masked_single_compare(&reconstructed_seq, &quals_per_base, alt_bytes, min_baseq);

        debug!(
            "Case B (ALT-only): reliable={} mismatches={}",
            reliable_count, mismatches
        );

        if reliable_count > 0 && mismatches == 0 {
            debug!("Matches ALT on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return (false, true, med_haplotype_qual);
        }
    } else if matches_ref_len {
        // Case C: Only REF length matches — no ambiguity possible
        let (mismatches, reliable_count) =
            masked_single_compare(&reconstructed_seq, &quals_per_base, ref_bytes, min_baseq);

        debug!(
            "Case C (REF-only): reliable={} mismatches={}",
            reliable_count, mismatches
        );

        if reliable_count > 0 && mismatches == 0 {
            debug!("Matches REF on {} reliable bases, med_qual={}", reliable_count, med_haplotype_qual);
            return (true, false, med_haplotype_qual);
        }
    } else {
        debug!(
            "Length mismatch: recon={} alt={} ref={}",
            recon_len,
            alt_bytes.len(),
            ref_bytes.len()
        );
    }

    // --- Phase 3: Alignment-based fallback (indelpost approach) ---
    // When narrow-window reconstruction fails (FM1: D-truncation, FM2: adjacent I),
    // expand to the full ref_context window and use dual-haplotype SW alignment.
    //
    // CRITICAL: Use raw read window extraction (not CIGAR-projected) to preserve
    // the true biological sequence. For complex variants (e.g. EPHA7 REF=TCC ALT=CT),
    // BWA represents ALT reads as DEL+INS CIGARs. CIGAR-projected extraction
    // produces a hybrid sequence matching neither REF nor ALT haplotype.
    // Raw extraction gives the contiguous read bases that SW can correctly align.
    //
    // Pre-filter (indelpost pattern): only attempt SW for reads showing
    // evidence of carrying the variant (soft-clips, indels near window).
    // This eliminates ~80-90% of clean REF reads from expensive alignment.
    if let Some(ref ctx) = variant.ref_context {
        let win_start = variant.ref_context_start;
        let win_end = win_start + ctx.len() as i64;

        if !is_worth_realignment(record, win_start, win_end) {
            trace!(
                "Phase 3 skipped: read has clean CIGAR over [{}, {})",
                win_start, win_end
            );
            return (false, false, 0);
        }

        let (sub_seq, sub_quals) = extract_raw_read_window(
            record, win_start, win_end, variant.pos, variant.ref_allele.len()
        );
        if sub_seq.len() >= 3 {
            debug!(
                "Phase 3 fallback: extracted {} raw bases over [{}, {})",
                sub_seq.len(), win_start, win_end
            );
            return classify_by_alignment(
                &sub_seq, &sub_quals, variant, min_baseq,
                alt_aligner, ref_aligner,
            );
        }
    }

    (false, false, 0)
}

/// Masked comparison against two alleles simultaneously.
///
/// Masks out bases with quality below `min_baseq` and counts mismatches against
/// both `allele_a` and `allele_b` on the remaining reliable bases only.
///
/// Returns `(mismatches_a, mismatches_b, reliable_count)`.
fn masked_dual_compare(
    recon: &[u8],
    quals: &[u8],
    allele_a: &[u8],
    allele_b: &[u8],
    min_baseq: u8,
) -> (usize, usize, usize) {
    let mut mm_a = 0;
    let mut mm_b = 0;
    let mut reliable = 0;

    for (i, &base) in recon.iter().enumerate() {
        if quals[i] < min_baseq {
            continue; // mask out low-quality base entirely
        }
        reliable += 1;
        let b = base.to_ascii_uppercase();
        if b != allele_a[i].to_ascii_uppercase() {
            mm_a += 1;
        }
        if b != allele_b[i].to_ascii_uppercase() {
            mm_b += 1;
        }
    }

    (mm_a, mm_b, reliable)
}

/// Masked comparison against a single allele.
///
/// Masks out bases below `min_baseq` and counts mismatches on reliable bases.
/// Any mismatch on a reliable base means no match.
///
/// Returns `(mismatches, reliable_count)`.
fn masked_single_compare(
    recon: &[u8],
    quals: &[u8],
    allele: &[u8],
    min_baseq: u8,
) -> (usize, usize) {
    let mut mismatches = 0;
    let mut reliable = 0;

    for (i, &base) in recon.iter().enumerate() {
        if quals[i] < min_baseq {
            continue;
        }
        reliable += 1;
        if !base.eq_ignore_ascii_case(&allele[i]) {
            mismatches += 1;
        }
    }

    (mismatches, reliable)
}

/// Compute the median quality of bases that pass the minimum threshold.
///
/// Follows the GATK `BaseQuality` annotation standard (median rather than min)
/// to prevent a single low-quality outlier from penalizing an entire read's
/// contribution to fragment consensus.
///
/// Returns 0 if no qualifying bases.
#[inline]
fn median_qual(quals: &[u8], min_baseq: u8) -> u8 {
    let mut filtered: Vec<u8> = quals.iter()
        .copied()
        .filter(|&q| q >= min_baseq)
        .collect();
    if filtered.is_empty() { return 0; }
    filtered.sort_unstable();
    filtered[filtered.len() / 2]
}

/// Returns (is_ref, is_alt, base_qual) for SNP variants.
/// Quality is the base quality at the variant position.
fn check_snp(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    let read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return (false, false, 0),
    };

    let qual = record.qual()[read_pos];
    if qual < min_baseq {
        return (false, false, 0);
    }

    let base = record.seq()[read_pos] as char;
    let ref_char = variant
        .ref_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let alt_char = variant
        .alt_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let base_upper = base.to_ascii_uppercase();

    let is_ref = base_upper == ref_char;
    let is_alt = base_upper == alt_char;
    // Return quality for fragment consensus scoring
    let base_qual = if is_ref || is_alt { qual } else { 0 };
    (is_ref, is_alt, base_qual)
}

/// Check if a read supports an insertion variant.
///
/// Returns (is_ref, is_alt, base_qual) where base_qual is the quality of the
/// anchor base, used for fragment-level consensus scoring.
///
/// Uses a single CIGAR walk with three detection strategies:
/// 1. **Backward boundary check:** When anchor falls at the start of an M block
///    and the previous CIGAR op was a matching insertion (fixes off-by-one at
///    M/I/M boundaries where the aligner splits the match block).
/// 2. **Strict match (fast path):** Insertion immediately after the anchor base.
///    Returns ALT immediately if length + sequence match.
/// 3. **Windowed scan (fallback):** Any insertion within ±5bp of the anchor,
///    validated by three safeguards:
///    - S1: Inserted sequence matches expected ALT bases (quality-masked)
///    - S2: Closest match wins (minimum |shift_pos - anchor_pos|)
///    - S3: Reference base at shifted anchor matches original anchor base
///      (via variant.ref_context)
/// 4. **Phase 3 haplotype fallback:** When a length-matching insertion exists
///    nearby but fails the sequence check (e.g., same biological event
///    represented differently by caller vs aligner), falls back to
///    check_complex for Smith-Waterman haplotype comparison.
fn check_insertion<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> (bool, bool, u8) {
    let cigar_view = record.cigar();
    let quals = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;

    let anchor_pos = variant.pos;
    let expected_ins_len = variant.alt_allele.len() - 1; // VCF ALT includes anchor
    let expected_ins_seq = &variant.alt_allele.as_bytes()[1..]; // ALT without anchor
    let original_anchor_base = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();

    // Windowed scan parameters
    let window: i64 = 5;
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    // State tracked across the CIGAR walk
    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None; // distance of best windowed match
    let mut has_nearby_length_match = false; // length-matching Ins found but seq check failed

    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let block_end = ref_pos + len_i64;

                // Track anchor read position whenever we encounter it
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    let offset = (anchor_pos - ref_pos) as usize;
                    anchor_read_pos = Some(read_pos + offset);
                }

                // --- P0-2: Backward boundary check ---
                // When anchor falls at block_end of prior M block, CIGAR geometry
                // places it at ref_pos of THIS block (after the Ins was consumed).
                // Check backward: was the previous op a matching insertion?
                if anchor_pos == ref_pos && i > 0 {
                    if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i - 1) {
                        let ins_len_usize = *ins_len as usize;
                        if ins_len_usize == expected_ins_len {
                            let ins_read_start = read_pos - ins_len_usize;
                            if ins_read_start + ins_len_usize <= record.seq().len() {
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_read_start..ins_read_start + ins_len_usize];
                                // P1-1: Quality-aware fuzzy match
                                let ins_quals = &quals[ins_read_start..ins_read_start + ins_len_usize];
                                let (mismatches, reliable) = masked_single_compare(
                                    ins_seq, ins_quals, expected_ins_seq, min_baseq
                                );
                                if reliable > 0 && mismatches == 0 {
                                    let qual = if read_pos < quals.len() { quals[read_pos] } else { 0 };
                                    debug!(
                                        "check_insertion: backward boundary match at pos {}, qual={}",
                                        anchor_pos, qual
                                    );
                                    return (false, true, qual); // ALT — backward match
                                }
                            }
                        }
                    }
                }

                // --- Strict fast path: anchor at end of this block ---
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    if anchor_pos == block_end - 1 {
                        // Anchor is the last base of this match block.
                        // Check if next op is an insertion with matching length + sequence.
                        if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                            let ins_len_usize = *ins_len as usize;
                            if ins_len_usize == expected_ins_len {
                                let ins_start = read_pos + *len as usize;
                                if ins_start + ins_len_usize <= record.seq().len() {
                                    let ins_seq = &record.seq().as_bytes()
                                        [ins_start..ins_start + ins_len_usize];
                                    // P1-1: Quality-aware fuzzy match
                                    let ins_quals = &quals[ins_start..ins_start + ins_len_usize];
                                    let (mismatches, reliable) = masked_single_compare(
                                        ins_seq, ins_quals, expected_ins_seq, min_baseq
                                    );
                                    if reliable > 0 && mismatches == 0 {
                                        let arp = anchor_read_pos.unwrap_or(0);
                                        let qual = if arp < quals.len() { quals[arp] } else { 0 };
                                        debug!(
                                            "check_insertion: strict match at pos {}, anchor_qual={}",
                                            anchor_pos, qual
                                        );
                                        return (false, true, qual); // ALT — strict match
                                    }
                                }
                            }
                        }
                        // Anchor at end but no matching insertion → REF coverage
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match block → read covers anchor without insertion
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check if any Ins after this block is within window ---
                if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                    let ins_ref_pos = block_end; // genomic position where insertion occurs
                    if ins_ref_pos >= window_start && ins_ref_pos <= window_end
                        && ins_ref_pos != anchor_pos + 1 // skip strict position (already handled)
                    {
                        let ins_len_usize = *ins_len as usize;
                        // Safeguard 1: length must match
                        if ins_len_usize == expected_ins_len {
                            let ins_start = read_pos + *len as usize;
                            if ins_start + ins_len_usize <= record.seq().len() {
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_start..ins_start + ins_len_usize];
                                // Quality-aware fuzzy match for inserted bases
                                let ins_quals = &quals[ins_start..ins_start + ins_len_usize];
                                let (mismatches, reliable) = masked_single_compare(
                                    ins_seq, ins_quals, expected_ins_seq, min_baseq
                                );
                                if reliable > 0 && mismatches == 0 {
                                    // Safeguard 3: verify anchor base at shifted position
                                    let shifted_anchor_pos = ins_ref_pos - 1;
                                    let anchor_ok = match &variant.ref_context {
                                        Some(ctx) => {
                                            let ctx_offset = (shifted_anchor_pos
                                                - variant.ref_context_start)
                                                as usize;
                                            if ctx_offset < ctx.len() {
                                                ctx.as_bytes()[ctx_offset].to_ascii_uppercase()
                                                    == original_anchor_base
                                            } else {
                                                debug!(
                                                    "ref_context offset {} out of bounds (len={}), allowing",
                                                    ctx_offset, ctx.len()
                                                );
                                                true
                                            }
                                        }
                                        None => true,
                                    };

                                    if anchor_ok {
                                        // Safeguard 2: track closest match
                                        let distance =
                                            (ins_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                        if best_windowed_match
                                            .is_none_or(|prev| distance < prev)
                                        {
                                            best_windowed_match = Some(distance);
                                        }
                                    } else {
                                        debug!(
                                            "check_insertion: S3 reject at shifted pos {} \
                                             (anchor base mismatch)",
                                            shifted_anchor_pos
                                        );
                                    }
                                } else {
                                    // Length matches but sequence differs — the caller
                                    // and aligner may represent the same event
                                    // differently (e.g., shifted insertion in a repeat).
                                    // Track this so Phase 3 SW can arbitrate.
                                    has_nearby_length_match = true;
                                    debug!(
                                        "check_insertion: windowed I({}) at pos {} seq \
                                         mismatch (mismatches={}, reliable={}), \
                                         flagging for Phase 3 fallback",
                                        ins_len_usize, ins_ref_pos,
                                        mismatches, reliable
                                    );
                                }
                            }
                        }
                    }
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            _ => {}
        }
    }

    // Anchor quality for fragment consensus (used for both ALT and REF returns)
    let anchor_qual = anchor_read_pos
        .filter(|&p| p < quals.len())
        .map(|p| quals[p])
        .unwrap_or(0);

    // Evaluate results after full CIGAR walk
    if best_windowed_match.is_some() {
        debug!(
            "check_insertion: windowed match for variant at pos {}, anchor_qual={}",
            anchor_pos, anchor_qual
        );
        return (false, true, anchor_qual); // ALT — windowed match
    }

    // Phase 3 haplotype fallback: when a length-matching insertion exists nearby
    // but the sequence check failed, the caller and aligner may represent the
    // same biological event differently (e.g., FLT4 A→GGAT placed at different
    // positions with different inserted bases). Suppress found_ref_coverage to
    // let check_complex → Smith-Waterman arbitrate using full haplotype comparison.
    if has_nearby_length_match && found_ref_coverage {
        debug!(
            "check_insertion: nearby I({}) with seq mismatch at pos {}, \
             falling back to check_complex for Phase 3 SW",
            expected_ins_len, anchor_pos
        );
        return check_complex(record, variant, min_baseq, alt_aligner, ref_aligner);
    }

    if found_ref_coverage {
        return (true, false, anchor_qual); // REF — read covers anchor without matching insertion
    }
    (false, false, 0) // Read does not cover the variant region
}

/// Check if a read supports a deletion variant.
///
/// Returns (is_ref, is_alt, base_qual) where base_qual is the quality of the
/// anchor base, used for fragment-level consensus scoring.
///
/// Uses the same single-walk strategy as check_insertion:
/// 1. **Strict match (fast path):** Deletion immediately after anchor, length matches.
/// 2. **Windowed scan (fallback):** Any deletion within ±5bp, validated by:
///    - S1: Deletion length matches expected
///    - S2: Closest match wins
///    - S3: Reference bases at shifted position match expected deleted sequence
///      (via variant.ref_context)
/// 3. **Haplotype fallback:** When CIGAR geometry doesn't match (e.g. different
///    breakpoint placement or wrong deletion length), delegates to `check_complex`
///    for quality-aware haplotype comparison.
fn check_deletion<F: Fn(u8, u8) -> i32>(
    record: &Record,
    variant: &Variant,
    min_baseq: u8,
    alt_aligner: &mut Aligner<F>,
    ref_aligner: &mut Aligner<F>,
) -> (bool, bool, u8) {
    let cigar_view = record.cigar();
    let quals = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos: usize = 0;

    let anchor_pos = variant.pos;
    let expected_del_len = variant.ref_allele.len() - 1; // REF without anchor
    // The expected deleted bases (REF without the anchor base)
    let expected_del_seq = &variant.ref_allele.as_bytes()[1..];

    // Windowed scan parameters
    let window: i64 = 5;
    let window_start = (anchor_pos - window).max(0);
    let window_end = anchor_pos + window;

    let mut found_ref_coverage = false;
    let mut anchor_read_pos: Option<usize> = None; // read position of anchor base
    let mut best_windowed_match: Option<u64> = None;
    let mut has_large_cigar_del = false; // tracks if read carries a large deletion

    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let block_end = ref_pos + len_i64;

                // Track anchor read position whenever we encounter it
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    let offset = (anchor_pos - ref_pos) as usize;
                    anchor_read_pos = Some(read_pos + offset);
                }

                // --- Strict fast path ---
                if anchor_pos >= ref_pos && anchor_pos < block_end {
                    if anchor_pos == block_end - 1 {
                        // Anchor at end of match. Check if next op is a deletion.
                        if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                            if *del_len as usize == expected_del_len {
                                let arp = anchor_read_pos.unwrap_or(0);
                                let qual = if arp < quals.len() { quals[arp] } else { 0 };
                                debug!(
                                    "check_deletion: strict match at pos {}, anchor_qual={}",
                                    anchor_pos, qual
                                );
                                return (false, true, qual); // ALT — strict match
                            } else {
                                // P0-3: D found at anchor but wrong length.
                                // Use SV-caller-style reciprocal overlap matching:
                                // aligners often report slightly different breakpoints
                                // for the same large deletion, producing different
                                // CIGAR D lengths. If both start at the same anchor
                                // and share ≥50% reciprocal overlap, treat as the
                                // same biological event.
                                // Precedent: SURVIVOR uses ≥50% overlap, BEDTools
                                // uses configurable reciprocal overlap for SV matching.
                                let found_del_len = *del_len as usize;
                                let min_del = expected_del_len.min(found_del_len);
                                let max_del = expected_del_len.max(found_del_len);
                                let reciprocal_overlap =
                                    min_del as f64 / max_del as f64;

                                if expected_del_len >= 50 && reciprocal_overlap >= 0.5 {
                                    // Large deletion with significant overlap → ALT
                                    let arp = anchor_read_pos.unwrap_or(0);
                                    let qual = if arp < quals.len() {
                                        quals[arp]
                                    } else {
                                        0
                                    };
                                    debug!(
                                        "check_deletion: tolerant match D({}) ≈ D({}) \
                                         at pos {} (reciprocal_overlap={:.2}, \
                                         threshold=0.50, anchor_qual={})",
                                        found_del_len,
                                        expected_del_len,
                                        anchor_pos,
                                        reciprocal_overlap,
                                        qual
                                    );
                                    return (false, true, qual); // ALT — tolerant match
                                }

                                // Small deletion or low overlap: fall back to
                                // check_complex for haplotype-based comparison.
                                debug!(
                                    "check_deletion: D({}) at anchor {} but expected D({}), \
                                     reciprocal_overlap={:.2} (below 0.50 or del<50bp), \
                                     falling back to check_complex",
                                    found_del_len,
                                    anchor_pos,
                                    expected_del_len,
                                    reciprocal_overlap
                                );
                                return check_complex(record, variant, min_baseq, alt_aligner, ref_aligner);
                            }
                        }
                        found_ref_coverage = true;
                    } else {
                        // Anchor in middle of match → REF coverage
                        found_ref_coverage = true;
                    }
                }

                // --- Windowed scan: check for Del after this block within window ---
                if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                    let del_ref_pos = block_end; // genomic position where deletion starts
                    if del_ref_pos >= window_start && del_ref_pos <= window_end
                        && del_ref_pos != anchor_pos + 1 // skip strict position
                    {
                        let del_len_usize = *del_len as usize;

                        // Safeguard 1: deletion length check.
                        // Accept exact matches, OR for large deletions (≥50bp),
                        // accept reciprocal overlap ≥50% (same logic as Fix 1
                        // on the strict path). Aligners often report slightly
                        // different breakpoints for the same biological event.
                        let length_ok = if del_len_usize == expected_del_len {
                            true
                        } else if expected_del_len >= 50 {
                            let min_del = expected_del_len.min(del_len_usize);
                            let max_del = expected_del_len.max(del_len_usize);
                            let overlap = min_del as f64 / max_del as f64;
                            if overlap >= 0.5 {
                                debug!(
                                    "check_deletion: windowed tolerant match \
                                     D({}) ≈ D({}) at pos {} (overlap={:.2})",
                                    del_len_usize, expected_del_len,
                                    del_ref_pos, overlap
                                );
                                // Track that this read carries a large deletion
                                has_large_cigar_del = true;
                                true
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        if length_ok {
                            // Safeguard 3: verify the deleted reference bases match
                            // (only for exact-length matches; skip for tolerant
                            // matches where lengths differ)
                            let del_ok = if del_len_usize != expected_del_len {
                                true // tolerant match — length differs, skip seq check
                            } else {
                                match &variant.ref_context {
                                    Some(ctx) => {
                                        let ctx_bytes = ctx.as_bytes();
                                        let ctx_offset =
                                            (del_ref_pos - variant.ref_context_start) as usize;
                                        if ctx_offset + del_len_usize <= ctx_bytes.len() {
                                            let ref_at_shift = &ctx_bytes
                                                [ctx_offset..ctx_offset + del_len_usize];
                                            let ref_quals = vec![u8::MAX; del_len_usize];
                                            let (mismatches, reliable) = masked_single_compare(
                                                ref_at_shift, &ref_quals, expected_del_seq, 0
                                            );
                                            reliable > 0 && mismatches == 0
                                        } else {
                                            debug!(
                                                "ref_context offset {} out of bounds (len={}), allowing",
                                                ctx_offset, ctx_bytes.len()
                                            );
                                            true
                                        }
                                    }
                                    None => true,
                                }
                            };

                            if del_ok {
                                // Safeguard 2: track closest match
                                let distance =
                                    (del_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                if best_windowed_match
                                    .is_none_or(|prev| distance < prev)
                                {
                                    best_windowed_match = Some(distance);
                                }
                            } else {
                                debug!(
                                    "check_deletion: S3 reject at shifted pos {} \
                                     (deleted bases mismatch)",
                                    del_ref_pos
                                );
                            }
                        }
                    }
                }

                ref_pos = block_end;
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            _ => {}
        }
    }

    // Anchor quality for fragment consensus (used for both ALT and REF returns)
    let anchor_qual = anchor_read_pos
        .filter(|&p| p < quals.len())
        .map(|p| quals[p])
        .unwrap_or(0);

    // Evaluate results after full CIGAR walk
    if best_windowed_match.is_some() {
        debug!(
            "check_deletion: windowed match for variant at pos {}, anchor_qual={}",
            anchor_pos, anchor_qual
        );
        return (false, true, anchor_qual); // ALT — windowed match
    }

    // Fix 3: Interior REF guard for large deletions.
    //
    // For large deletions (≥50bp), reads that start *inside* the expected
    // deletion span are reference reads — they align to sequence that would
    // be absent in the ALT allele. These reads won't cover the anchor
    // (found_ref_coverage == false) and have no matching deletion CIGAR,
    // so without this guard they fall back to check_complex → Phase 3
    // Smith-Waterman alignment.
    //
    // In Phase 3, the REF haplotype (e.g. 1044bp for a 1023bp deletion)
    // far exceeds the read length (~100bp), causing semiglobal alignment
    // to penalize it with gaps. The short ALT haplotype (~21bp) fits
    // easily, producing a higher score → false ALT call.
    //
    // Fix: detect reads starting inside the deletion span and classify
    // them as REF directly. Only applied for large deletions (≥50bp)
    // where the semiglobal alignment length mismatch is significant.
    //
    // Quality note: these interior reads don't cover the anchor base,
    // so anchor_qual is 0. Using 0 would cause FragmentEvidence::resolve()
    // to discard this fragment (has_ref = best_ref_qual > 0 fails),
    // silently inflating ALT VAF. Instead, use the read's median base
    // quality as proxy evidence — the read's physical presence inside the
    // deletion span IS the evidence for REF support.
    let read_start = record.pos();
    let del_region_end = anchor_pos + expected_del_len as i64;

    if expected_del_len >= 50
        && read_start > anchor_pos
        && read_start < del_region_end
        && !found_ref_coverage
        && !has_large_cigar_del  // exclude reads carrying the actual deletion
    {
        // Use median read quality as proxy since anchor isn't covered.
        // Floor at 1 to guarantee FragmentEvidence counts this fragment.
        let proxy_qual = median_qual(record.qual(), min_baseq);
        let interior_qual = if proxy_qual > 0 { proxy_qual } else { 1 };
        debug!(
            "check_deletion: read at pos {} starts inside deletion span \
             [{}, {}), calling REF (interior read, del_len={}, proxy_qual={})",
            read_start,
            anchor_pos + 1,
            del_region_end,
            expected_del_len,
            interior_qual
        );
        return (true, false, interior_qual);
    }

    // P0-3: Haplotype fallback — when strict/windowed CIGAR matching found no
    // deletion match and the read doesn't cover the anchor, try check_complex
    // which reconstructs the read's haplotype and does quality-aware comparison.
    // Only fall back when NOT found_ref_coverage to avoid false positives on
    // reads that genuinely show REF at this position.
    if !found_ref_coverage && best_windowed_match.is_none() {
        debug!(
            "check_deletion: no CIGAR match at pos {}, falling back to check_complex",
            anchor_pos
        );
        return check_complex(record, variant, min_baseq, alt_aligner, ref_aligner);
    }

    if found_ref_coverage {
        return (true, false, anchor_qual); // REF
    }
    (false, false, 0) // Read does not cover the variant region
}

/// Returns (is_ref, is_alt, base_qual) for MNP variants.
/// Quality is the median base quality across all positions in the MNP.
fn check_mnp(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    // MNP: REF=AT, ALT=CG. Lengths equal, > 1.
    let len = variant.ref_allele.len();
    
    // Find read position of the first base
    let start_read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return (false, false, 0),
    };

    // Check if the read covers the entire MNP
    if start_read_pos + len > record.seq().len() {
        return (false, false, 0);
    }

    // Check qualities and sequence
    let quals = record.qual();
    let seq = record.seq();
    let seq_bytes = seq.as_bytes();
    
    let ref_bytes = variant.ref_allele.as_bytes();
    let alt_bytes = variant.alt_allele.as_bytes();

    let mut matches_ref = true;
    let mut matches_alt = true;
    let mut mnp_quals: Vec<u8> = Vec::with_capacity(len);

    for i in 0..len {
        let pos = start_read_pos + i;
        
        // Check quality
        if quals[pos] < min_baseq {
            return (false, false, 0);
        }
        mnp_quals.push(quals[pos]);

        let base = seq_bytes[pos].to_ascii_uppercase();
        let r = ref_bytes[i].to_ascii_uppercase();
        let a = alt_bytes[i].to_ascii_uppercase();

        if base != r {
            matches_ref = false;
        }
        if base != a {
            matches_alt = false;
        }
    }

    // Ensure no indels in the MNP region
    let end_read_pos = match find_read_pos(record, variant.pos + len as i64 - 1) {
        Some(p) => p,
        None => return (false, false, 0),
    };

    // If contiguous, end - start should be len - 1
    if end_read_pos - start_read_pos != len - 1 {
        return (false, false, 0); // Indel detected within MNP
    }

    // Return quality only for matching alleles
    let med_qual = if matches_ref || matches_alt { median_qual(&mnp_quals, min_baseq) } else { 0 };
    (matches_ref, matches_alt, med_qual)
}

/// Find the read index corresponding to a genomic position.
fn find_read_pos(record: &Record, target_pos: i64) -> Option<usize> {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return Some(read_pos + (target_pos - ref_pos) as usize);
                }
                ref_pos += *len as i64;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return None; // Position is deleted
                }
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}
