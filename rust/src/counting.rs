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
use log::debug;

/// Evidence accumulated for a single fragment (read pair) at a variant site.
/// Tracks the best base quality seen for each allele across both reads,
/// enabling quality-weighted consensus to resolve R1-vs-R2 conflicts.
#[derive(Debug, Clone)]
struct FragmentEvidence {
    /// Best base quality seen supporting REF across reads in this fragment
    best_ref_qual: u8,
    /// Best base quality seen supporting ALT across reads in this fragment
    best_alt_qual: u8,
    /// Orientation of read 1 (true = forward, false = reverse), if seen
    read1_orientation: Option<bool>,
    /// Orientation of read 2 (true = forward, false = reverse), if seen
    read2_orientation: Option<bool>,
}

impl FragmentEvidence {
    fn new() -> Self {
        FragmentEvidence {
            best_ref_qual: 0,
            best_alt_qual: 0,
            read1_orientation: None,
            read2_orientation: None,
        }
    }

    /// Record a read's allele call and orientation into this fragment's evidence.
    fn observe(&mut self, is_ref: bool, is_alt: bool, base_qual: u8, is_read1: bool, is_forward: bool) {
        if is_ref && base_qual > self.best_ref_qual {
            self.best_ref_qual = base_qual;
        }
        if is_alt && base_qual > self.best_alt_qual {
            self.best_alt_qual = base_qual;
        }
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

    /// Get the fragment orientation. Prefers read 1 orientation when both reads present.
    fn orientation(&self) -> Option<bool> {
        // Prefer read 1 orientation (consistent with original behavior)
        self.read1_orientation.or(self.read2_orientation)
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
#[pyfunction]
#[pyo3(signature = (bam_path, variants, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10))]
pub fn count_bam(
    py: Python<'_>,
    bam_path: String,
    variants: Vec<Variant>,
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

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<BaseCounts>, anyhow::Error> = py.allow_threads(move || {
        pool.install(|| {
            variants
                .par_iter()
                .map_init(
                    || {
                        // Initialize thread-local BAM reader
                        bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                            anyhow::anyhow!("Failed to open BAM: {}", e)
                        })
                    },
                    |bam_result, variant| {
                        // Get the reader or return error if initialization failed
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        count_single_variant(
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
                        )
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
        let (is_ref, is_alt, base_qual) = check_allele_with_qual(&record, variant, min_baseq);

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
    for (_qname_hash, evidence) in &fragments {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);

        // Get fragment orientation (prefer read 1)
        let orientation = match evidence.orientation() {
            Some(o) => o,
            None => continue, // No orientation data (should not happen)
        };

        // Count every fragment in dpf regardless of consensus outcome.
        // Discarded fragments (ambiguous R1-vs-R2 within quality threshold)
        // are still real molecules — tracking them in dpf makes the gap
        // dpf - (rdf + adf) a useful quality metric for the locus.
        counts.dpf += 1;

        if frag_ref {
            counts.rdf += 1;
            if orientation {
                counts.rdf_fwd += 1;
            } else {
                counts.rdf_rev += 1;
            }
        } else if frag_alt {
            counts.adf += 1;
            if orientation {
                counts.adf_fwd += 1;
            } else {
                counts.adf_rev += 1;
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
fn check_allele_with_qual(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    let variant_type = &variant.variant_type;
    debug!("check_allele type={} pos={} ref={} alt={}", variant_type, variant.pos, variant.ref_allele, variant.alt_allele);

    match variant_type.as_str() {
        "SNP" => check_snp(record, variant, min_baseq),
        "INSERTION" => check_insertion(record, variant, min_baseq),
        "DELETION" => check_deletion(record, variant, min_baseq),
        "MNP" => {
            debug!("Dispatching to check_mnp");
            check_mnp(record, variant, min_baseq)
        },
        "COMPLEX" => check_complex(record, variant, min_baseq),
        _ => {
            // Auto-detect MNP if type is not explicit but looks like one
            if variant.ref_allele.len() == variant.alt_allele.len()
                && variant.ref_allele.len() > 1
            {
                debug!("Auto-detected MNP");
                check_mnp(record, variant, min_baseq)
            } else {
                // Fallback to complex check for anything else (e.g. DelIns)
                check_complex(record, variant, min_baseq)
            }
        }
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
/// Returns (is_ref, is_alt, base_qual) where base_qual is the minimum quality
/// across the reconstructed haplotype bases, used for fragment consensus.
fn check_complex(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
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
                read_pos += *len as usize;
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

    // Minimum quality across the reconstructed haplotype for fragment consensus.
    // This is the quality we return for whichever allele matches.
    let min_haplotype_qual = quals_per_base.iter().copied().min().unwrap_or(0);

    // --- Phase 2: Quality-Aware Masked Comparison ---
    // Mask out low-quality bases. Only reliable bases (qual >= min_baseq) vote.
    let alt_bytes = variant.alt_allele.as_bytes();
    let ref_bytes = variant.ref_allele.as_bytes();
    let recon_len = reconstructed_seq.len();
    let matches_alt_len = recon_len == alt_bytes.len();
    let matches_ref_len = recon_len == ref_bytes.len();

    if matches_alt_len && matches_ref_len {
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
            debug!("Matches ALT on {} reliable bases, min_qual={}", reliable_count, min_haplotype_qual);
            return (false, true, min_haplotype_qual);
        }
        if mismatches_ref == 0 {
            debug!("Matches REF on {} reliable bases, min_qual={}", reliable_count, min_haplotype_qual);
            return (true, false, min_haplotype_qual);
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
            debug!("Matches ALT on {} reliable bases, min_qual={}", reliable_count, min_haplotype_qual);
            return (false, true, min_haplotype_qual);
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
            debug!("Matches REF on {} reliable bases, min_qual={}", reliable_count, min_haplotype_qual);
            return (true, false, min_haplotype_qual);
        }
    } else {
        debug!(
            "Length mismatch: recon={} alt={} ref={}",
            recon_len,
            alt_bytes.len(),
            ref_bytes.len()
        );
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
        if base.to_ascii_uppercase() != allele[i].to_ascii_uppercase() {
            mismatches += 1;
        }
    }

    (mismatches, reliable)
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
/// Uses a single CIGAR walk with two detection strategies:
/// 1. **Strict match (fast path):** Insertion immediately after the anchor base.
///    Returns ALT immediately if length + sequence match.
/// 2. **Windowed scan (fallback):** Any insertion within ±5bp of the anchor,
///    validated by three safeguards:
///    - S1: Inserted sequence matches expected ALT bases exactly
///    - S2: Closest match wins (minimum |shift_pos - anchor_pos|)
///    - S3: Reference base at shifted anchor matches original anchor base
///          (via variant.ref_context)
fn check_insertion(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    let _ = min_baseq; // Quality is read for fragment scoring, not filtering (CIGAR-based detection)
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
                                    if ins_seq == expected_ins_seq {
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
                        // Safeguard 1: sequence identity — length and bases must match
                        if ins_len_usize == expected_ins_len {
                            let ins_start = read_pos + *len as usize;
                            if ins_start + ins_len_usize <= record.seq().len() {
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_start..ins_start + ins_len_usize];
                                if ins_seq == expected_ins_seq {
                                    // Safeguard 3: verify anchor base at shifted position
                                    // The base before the insertion (ins_ref_pos - 1) should
                                    // match the original anchor base
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
                                                false // offset out of bounds
                                            }
                                        }
                                        None => true, // no ref_context → skip check
                                    };

                                    if anchor_ok {
                                        // Safeguard 2: track closest match
                                        let distance =
                                            (ins_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                        if best_windowed_match
                                            .map_or(true, |prev| distance < prev)
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
///          (via variant.ref_context)
fn check_deletion(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    let _ = min_baseq; // Quality is read for fragment scoring, not filtering (CIGAR-based detection)
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
                        // Safeguard 1: deletion length matches
                        if del_len_usize == expected_del_len {
                            // Safeguard 3: verify the deleted reference bases match
                            let del_ok = match &variant.ref_context {
                                Some(ctx) => {
                                    let ctx_bytes = ctx.as_bytes();
                                    let ctx_offset =
                                        (del_ref_pos - variant.ref_context_start) as usize;
                                    if ctx_offset + del_len_usize <= ctx_bytes.len() {
                                        let ref_at_shift = &ctx_bytes
                                            [ctx_offset..ctx_offset + del_len_usize];
                                        // Compare case-insensitively
                                        ref_at_shift.iter().zip(expected_del_seq.iter()).all(
                                            |(a, b)| {
                                                a.to_ascii_uppercase()
                                                    == b.to_ascii_uppercase()
                                            },
                                        )
                                    } else {
                                        false // offset out of bounds
                                    }
                                }
                                None => true, // no ref_context → skip check
                            };

                            if del_ok {
                                // Safeguard 2: track closest match
                                let distance =
                                    (del_ref_pos - (anchor_pos + 1)).unsigned_abs();
                                if best_windowed_match
                                    .map_or(true, |prev| distance < prev)
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
    if found_ref_coverage {
        return (true, false, anchor_qual); // REF
    }
    (false, false, 0) // Read does not cover the variant region
}

/// Returns (is_ref, is_alt, base_qual) for MNP variants.
/// Quality is the minimum base quality across all positions in the MNP.
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
    let mut min_qual: u8 = u8::MAX; // Track minimum quality across MNP positions

    for i in 0..len {
        let pos = start_read_pos + i;
        
        // Check quality
        if quals[pos] < min_baseq {
            return (false, false, 0);
        }
        // Track minimum quality for fragment consensus
        if quals[pos] < min_qual {
            min_qual = quals[pos];
        }

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
    let base_qual = if matches_ref || matches_alt { min_qual } else { 0 };
    (matches_ref, matches_alt, base_qual)
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
