//! Base counting engine for variant classification.
//!
//! This module orchestrates read-level allele classification for a list of
//! variants in a BAM file. It coordinates the variant-type-specific checkers,
//! fragment-level evidence tracking, and alignment backends.
//!
//! ## Submodules
//!
//! - [`fragment`] — Fragment evidence tracking and QNAME hashing
//! - [`alignment`] — Smith-Waterman alignment backend (Phase 3)
//! - [`pairhmm`] — PairHMM alignment backend (probabilistic Phase 3 alternative)
//! - [`variant_checks`] — Per-variant-type classification (SNP, MNP, Ins, Del, Complex)
//! - [`utils`] — Shared utility functions (position lookup, quality, masked comparison)

mod fragment;
mod alignment;
pub mod pairhmm;
mod variant_checks;
mod utils;

use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Record};
use std::collections::HashMap;

use crate::stats::fisher_strand_bias;
use crate::types::{BaseCounts, Variant};

use rayon::prelude::*;

use anyhow::{Context, Result};
use log::debug;
use bio::alignment::pairwise::Aligner;

use fragment::{FragmentEvidence, hash_qname};
use variant_checks::{check_snp, check_mnp, check_complex, check_insertion, check_deletion, MnpResult};


/// Alignment backend selection for Phase 3 fallback classification.
///
/// Controls which algorithm is used when variant-type-specific checkers
/// (SNP, Ins, Del, MNP, Complex) need to fall back to haplotype-level
/// alignment for ambiguous reads.
///
/// Selectable via `--alignment-backend` CLI flag.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum AlignmentBackend {
    /// Smith-Waterman with affine gap penalties (default for v2.8.0).
    /// Score-margin classification: (alt_score - ref_score) > 0 → ALT.
    #[default]
    SmithWaterman,
    /// PairHMM with BQ-aware emissions (planned default for v3.0.0).
    /// LLR classification: log P(read|ALT) - log P(read|REF) > threshold → ALT.
    PairHMM {
        /// Log-likelihood ratio threshold for confident calls (default: 2.3 ≈ 10:1 odds).
        llr_threshold: f64,
        /// Gap-open probability (linear scale, default: 1e-4).
        gap_open: f64,
        /// Gap-extend probability (linear scale, default: 0.1).
        gap_extend: f64,
        /// Gap-open probability for repeat regions (linear scale, default: 1e-2).
        gap_open_repeat: f64,
        /// Gap-extend probability for repeat regions (linear scale, default: 0.5).
        gap_extend_repeat: f64,
    },
}

impl AlignmentBackend {
    /// Create PairHMM backend with default parameters.
    ///
    /// Convenience constructor for tests and downstream consumers.
    /// In production, the Python CLI passes params directly to the
    /// `PairHMM { ... }` variant via `count_bam()`.
    #[allow(dead_code)]
    pub fn pairhmm_default() -> Self {
        AlignmentBackend::PairHMM {
            llr_threshold: 2.3,
            gap_open: 1e-4,
            gap_extend: 0.1,
            gap_open_repeat: 1e-2,
            gap_extend_repeat: 0.5,
        }
    }
}


/// Count bases for a list of variants in a BAM file.
///
/// When `decomposed` is provided (same length as `variants`), variants with
/// a `Some(decomposed_variant)` are counted twice — once with the original
/// allele and once with the corrected allele. The result with the higher
/// `ad` (alt_count) is returned, with `used_decomposed` set accordingly.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (bam_path, variants, decomposed, min_mapq, min_baseq, filter_duplicates, filter_secondary, filter_supplementary, filter_qc_failed, filter_improper_pair, filter_indel, threads, fragment_qual_threshold=10, sibling_variants=Vec::new(), alignment_backend="sw", hmm_llr_threshold=2.3, hmm_gap_open=1e-4, hmm_gap_extend=0.1, hmm_gap_open_repeat=1e-2, hmm_gap_extend_repeat=0.5))]
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
    sibling_variants: Vec<Vec<Variant>>,
    alignment_backend: &str,
    hmm_llr_threshold: f64,
    hmm_gap_open: f64,
    hmm_gap_extend: f64,
    hmm_gap_open_repeat: f64,
    hmm_gap_extend_repeat: f64,
) -> PyResult<Vec<BaseCounts>> {
    // Parse alignment backend from string
    let backend = match alignment_backend {
        "hmm" | "pairhmm" => AlignmentBackend::PairHMM {
            llr_threshold: hmm_llr_threshold,
            gap_open: hmm_gap_open,
            gap_extend: hmm_gap_extend,
            gap_open_repeat: hmm_gap_open_repeat,
            gap_extend_repeat: hmm_gap_extend_repeat,
        },
        _ => AlignmentBackend::SmithWaterman,
    };

    // We cannot share a single IndexedReader across threads because it's not Sync.
    // Instead, we use rayon's map_init to initialize a reader for each thread.
    // This is efficient because map_init reuses the thread-local state (the reader)
    // for multiple items processed by that thread.

    // Configure thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to build thread pool: {}", e)))?;

    // Pad sibling_variants to match variants length (handles default empty case)
    let mut sibling_variants = sibling_variants;
    let n = variants.len();
    sibling_variants.resize_with(n, Vec::new);

    // Zip variants with their decomposed counterparts and sibling alts for parallel iteration
    let paired: Vec<_> = variants.into_iter()
        .zip(decomposed)
        .zip(sibling_variants)
        .map(|((v, d), s)| (v, d, s))
        .collect();

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
                    |bam_result, (variant, decomp_opt, siblings)| {
                        // Get the reader or return error if initialization failed
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        let counts_orig = count_single_variant(
                            bam,
                            variant,
                            siblings,
                            min_mapq,
                            min_baseq,
                            filter_duplicates,
                            filter_secondary,
                            filter_supplementary,
                            filter_qc_failed,
                            filter_improper_pair,
                            filter_indel,
                            fragment_qual_threshold,
                            &backend,
                        )?;

                        // Dual-count: if a decomposed variant exists, count it too
                        // and return whichever has the higher alt_count.
                        if let Some(decomp) = decomp_opt {
                            let counts_decomp = count_single_variant(
                                bam,
                                decomp,
                                siblings,
                                min_mapq,
                                min_baseq,
                                filter_duplicates,
                                filter_secondary,
                                filter_supplementary,
                                filter_qc_failed,
                                filter_improper_pair,
                                filter_indel,
                                fragment_qual_threshold,
                                &backend,
                            )?;

                            if counts_decomp.ad > counts_orig.ad {
                                // Sanity: both hypotheses count the same reads at the
                                // same locus, so DP should be nearly identical. A large
                                // divergence indicates a counting bug.
                                if (counts_decomp.dp as i64 - counts_orig.dp as i64).abs() > 2 {
                                    log::warn!(
                                        "DP mismatch in dual-counting: decomp={} orig={} at {}:{} {}→{}",
                                        counts_decomp.dp, counts_orig.dp,
                                        variant.chrom, variant.pos + 1,
                                        variant.ref_allele, decomp.alt_allele
                                    );
                                }
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
    sibling_variants: &[Variant],
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    fragment_qual_threshold: u8,
    backend: &AlignmentBackend,
) -> Result<BaseCounts> {
    let tid = bam.header().tid(variant.chrom.as_bytes()).ok_or_else(|| {
        anyhow::anyhow!("Chromosome not found in BAM: {}", variant.chrom)
    })?;

    // Fetch region around the variant. For windowed indel detection,
    // we expand the window so that reads with shifted indels are also retrieved.
    // The window scales with repeat_span to capture indels that aligners
    // shift beyond 5bp in long homopolymers/microsatellites.
    let window_pad: i64 = std::cmp::max(5, variant.repeat_span as i64 + 2);
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
    // Dynamic gap_extend: in tandem repeat regions (repeat_span >= 10bp),
    // gap extension is free (0) to absorb polymerase slippage noise.
    // For stable DNA, gap_extend = -1 (identical to previous behavior).
    let gap_open: i32 = -5;
    let gap_extend: i32 = if variant.repeat_span >= 10 { 0 } else { -1 };
    let mut alt_aligner = Aligner::new(gap_open, gap_extend, &score_fn);
    let mut ref_aligner = Aligner::new(gap_open, gap_extend, &score_fn);

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
            &record, variant, min_baseq, &mut alt_aligner, &mut ref_aligner, backend,
        );

        // ── ANCHOR OVERLAP CHECK: The BAM fetch region is wider than the
        // variant footprint (±5bp for shifted indel detection). Reads that
        // were brought in by the wider fetch but don't actually overlap the
        // variant's anchor position should NOT count toward DP — unless they
        // were classified as REF/ALT (e.g., shifted indel matched via
        // windowed detection). This ensures DP matches samtools pileup at
        // the variant position.
        let ref_len = variant.ref_allele.len().max(1) as i64;
        let read_start = record.pos();
        let read_end = {
            let mut rend = record.pos();
            for op in record.cigar().iter() {
                match op {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len)
                    | Cigar::Del(len) | Cigar::RefSkip(len) => {
                        rend += *len as i64;
                    }
                    _ => {}
                }
            }
            rend
        };
        let overlaps_anchor = read_start < (variant.pos + ref_len)
            && read_end > variant.pos;
        if !overlaps_anchor && !is_ref && !is_alt {
            continue;
        }

        // ── TOTAL DEPTH: count all reads that overlap the variant anchor,
        // regardless of allele classification. This ensures DP reflects true
        // physical coverage even for reads that are neither REF nor ALT
        // (e.g., duplex N bases, third alleles at multi-allelic sites).
        counts.dp += 1;
        let is_reverse = record.is_reverse();

        if is_reverse {
            counts.dp_rev += 1;
        } else {
            counts.dp_fwd += 1;
        }

        // ── FRAGMENT TRACKING: track ALL fragments for DPF.
        // FragmentEvidence::observe() correctly handles (false, false) —
        // it skips updating best_ref_qual/best_alt_qual but still tracks
        // the fragment for DPF in the downstream resolution loop.
        let qname_hash = hash_qname(record.qname());
        let is_read1 = record.is_first_in_template();
        let is_forward = !is_reverse;

        let evidence = fragments.entry(qname_hash).or_insert_with(FragmentEvidence::new);
        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward);

        // ── ALLELE-SPECIFIC COUNTS: only REF/ALT reads contribute to RD/AD.
        // DP and DPF are already recorded above.
        if !is_ref && !is_alt {
            continue;
        }

        if is_ref {
            // Multi-allelic guard: if this read is classified as ALT for any
            // sibling variant at this locus, don't count it as REF for this
            // variant. This handles overlapping indels/complex variants where
            // a read carrying one variant's ALT could be miscounted as REF
            // for another variant at the same locus.
            if !sibling_variants.is_empty() {
                let mut is_sibling_alt = false;
                for sib in sibling_variants {
                    let (_sib_ref, sib_alt, _sib_qual) = check_allele_with_qual(
                        &record, sib, min_baseq,
                        &mut alt_aligner, &mut ref_aligner, backend,
                    );
                    if sib_alt {
                        debug!(
                            "Multi-allelic guard: read is ALT for sibling {}>{} at {}:{}, \
                             excluding from REF for {}>{}",
                            sib.ref_allele, sib.alt_allele,
                            variant.chrom, variant.pos + 1,
                            variant.ref_allele, variant.alt_allele,
                        );
                        is_sibling_alt = true;
                        break;
                    }
                }
                if is_sibling_alt {
                    continue; // Skip REF counting — this read belongs to a sibling
                }
            }
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
    backend: &AlignmentBackend,
) -> (bool, bool, u8) {
    // Dispatch based on allele lengths rather than the variant_type string.
    // This is more robust than relying on upstream type labels, which can be
    // inconsistent (e.g., a caller emitting "COMPLEX" for what is really a
    // pure deletion after normalization).
    let ref_len = variant.ref_allele.len();
    let alt_len = variant.alt_allele.len();
    debug!(
        "check_allele ref_len={} alt_len={} pos={} ref={} alt={}",
        ref_len, alt_len, variant.pos, variant.ref_allele, variant.alt_allele
    );

    if ref_len == 1 && alt_len == 1 {
        // SNP: single base substitution — no Phase 3 needed
        check_snp(record, variant, min_baseq)
    } else if ref_len == alt_len {
        // MNP: min-BQ-across-block strategy with selective Phase 3 fallback.
        match check_mnp(record, variant, min_baseq) {
            MnpResult::Ref(q) => (true, false, q),
            MnpResult::Alt(q) => (false, true, q),
            MnpResult::LowQuality => (false, false, 0),
            MnpResult::ThirdAllele => (false, false, 0),
            MnpResult::Structural => {
                debug!("MNP structural issue, falling back to Phase 3");
                check_complex(record, variant, min_baseq, alt_aligner, ref_aligner, backend)
            }
        }
    } else if ref_len == 1 {
        // Pure insertion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
        check_insertion(record, variant, min_baseq, alt_aligner, ref_aligner, backend)
    } else if alt_len == 1 {
        // Pure deletion: CIGAR-based fast paths, then backend-aware Phase 3 fallback
        check_deletion(record, variant, min_baseq, alt_aligner, ref_aligner, backend)
    } else {
        // Complex: ref_len != alt_len, both > 1 (e.g., DelIns)
        check_complex(record, variant, min_baseq, alt_aligner, ref_aligner, backend)
    }
}




#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use std::ffi::CString;

    /// Build a synthetic BAM record for testing.
    ///
    /// Creates a minimal Record with the given sequence, qualities, CIGAR,
    /// and mapping position. All other fields are set to sensible defaults
    /// (unmapped=false, mapq=60, proper_pair=true, forward strand).
    fn build_record(seq: &[u8], qual: &[u8], cigar: &CigarString, pos: i64) -> Record {
        let mut record = Record::new();
        let qname = CString::new("test_read").unwrap();
        record.set(
            qname.as_bytes(), // qname
            Some(cigar),       // cigar
            seq,               // seq
            qual,              // qual
        );
        record.set_pos(pos);
        record.set_tid(0);
        record.set_mapq(60);
        // Set as mapped, proper pair, first in template
        record.set_flags(0x01 | 0x02 | 0x40); // paired + proper + read1
        record
    }

    /// Helper: build a Variant with the given position, ref, and alt alleles.
    fn build_variant(pos: i64, ref_allele: &str, alt_allele: &str) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: String::new(),
            ref_context: None,
            ref_context_start: 0,
            repeat_span: 0,
        }
    }

    // ── MNP check_mnp tests ──

    #[test]
    fn test_mnp_high_quality_matches_ref() {
        // 2bp MNP: REF=AT, ALT=CG. Read has AT (matches REF).
        // All qualities >= 20.
        let seq = b"AAATGCC";
        let qual = &[30, 30, 30, 35, 30, 30, 30]; // all high quality
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant at pos 102-103 (0-based), REF=AT, ALT=CG
        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        match result {
            MnpResult::Ref(q) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
            other => panic!("Expected MnpResult::Ref, got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_high_quality_matches_alt() {
        // 2bp MNP: REF=AT, ALT=CG. Read has CG (matches ALT).
        let seq = b"AACGGCC";
        let qual = &[30, 30, 35, 38, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        match result {
            MnpResult::Alt(q) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
            other => panic!("Expected MnpResult::Alt, got {:?}", format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_one_low_quality_base() {
        // 2bp MNP: REF=AT, ALT=CG. Read has CG (matches ALT).
        // But one base has Q=8 < min_baseq=20.
        // min(35, 8) = 8 < 20 → LowQuality, NOT Phase 3 fallback.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 35, 8, 30, 30, 30]; // pos 3: Q=8 < 20
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::LowQuality),
            "Expected LowQuality for min(BQ) < threshold"
        );
    }

    #[test]
    fn test_mnp_all_low_quality() {
        // Both bases below threshold.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 5, 8, 30, 30, 30]; // pos 2: Q=5, pos 3: Q=8
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::LowQuality),
            "Expected LowQuality when all bases below threshold"
        );
    }

    #[test]
    fn test_mnp_third_allele() {
        // 2bp MNP: REF=AT, ALT=CG. Read has TT (matches neither).
        // All qualities pass threshold.
        let seq = b"AATTGCC";
        let qual = &[30, 30, 35, 38, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::ThirdAllele),
            "Expected ThirdAllele when bases match neither REF nor ALT"
        );
    }

    #[test]
    fn test_mnp_indel_in_block_structural() {
        // 3bp MNP: REF=ATG, ALT=CGC. Read has CIGAR 1M 1I 1D 2M,
        // creating a non-contiguous mapping across the MNP block.
        // The contiguity check should detect this and return Structural.
        let seq = b"AACXGCC"; // 7 bases, but CIGAR rearranges alignment
        let qual = &[30, 30, 35, 38, 35, 30, 30];
        // CIGAR: 2M 1I 1D 3M = consumes 2+0+1+3=6 ref, 2+1+0+3=6 read
        let cigar = CigarString(vec![
            Cigar::Match(2),
            Cigar::Ins(1),
            Cigar::Del(1),
            Cigar::Match(3),
        ]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant spans pos 102-104 (3bp MNP)
        // With the indel, positions won't be contiguous in read space
        let variant = build_variant(102, "ATG", "CGC");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when indel exists in MNP block"
        );
    }

    #[test]
    fn test_mnp_partial_coverage_structural() {
        // 3bp MNP but read only covers first 2 positions.
        // Read: 3 bases starting at pos 101, so covers 101-103.
        // Variant: pos 102-104 → pos 104 not covered.
        let seq = b"ACG";
        let qual = &[30, 35, 38];
        let cigar = CigarString(vec![Cigar::Match(3)]);
        let record = build_record(seq, qual, &cigar, 101);

        let variant = build_variant(102, "ATG", "CGC");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when read doesn't fully cover MNP"
        );
    }

    #[test]
    fn test_mnp_position_not_found_structural() {
        // Read starts after the variant position.
        let seq = b"ATCGATCG";
        let qual = &[30; 8];
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let record = build_record(seq, qual, &cigar, 200); // far beyond variant

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        assert!(
            matches!(result, MnpResult::Structural),
            "Expected Structural when variant position not in read"
        );
    }

    #[test]
    fn test_mnp_5bp_tert_pattern() {
        // 5bp MNP mimicking TERT: GAGGG→AAGGA
        // Read carries the ALT allele AAGGA at high quality.
        let seq = b"CCAAGGATTT";
        let qual = &[30, 30, 35, 38, 32, 36, 34, 30, 30, 30];
        let cigar = CigarString(vec![Cigar::Match(10)]);
        let record = build_record(seq, qual, &cigar, 100);

        // Variant at pos 102-106
        let variant = build_variant(102, "GAGGG", "AAGGA");
        let result = check_mnp(&record, &variant, 20);

        match result {
            MnpResult::Alt(q) => assert!(q >= 20, "Quality should be >= 20, got {}", q),
            other => panic!("Expected MnpResult::Alt for TERT-like 5bp MNP, got {:?}",
                           format_mnp_result(&other)),
        }
    }

    #[test]
    fn test_mnp_min_bq_boundary() {
        // Test the exact boundary: min(BQ) == min_baseq should PASS.
        let seq = b"AACGGCC";
        let qual = &[30, 30, 20, 30, 30, 30, 30]; // min BQ = 20 == threshold
        let cigar = CigarString(vec![Cigar::Match(7)]);
        let record = build_record(seq, qual, &cigar, 100);

        let variant = build_variant(102, "AT", "CG");
        let result = check_mnp(&record, &variant, 20);

        match result {
            MnpResult::Alt(q) => assert!(q >= 20, "Should pass at exact threshold boundary"),
            other => panic!("Expected Alt at exact BQ threshold, got {:?}",
                           format_mnp_result(&other)),
        }
    }

    /// Helper to format MnpResult for panic messages
    fn format_mnp_result(result: &MnpResult) -> String {
        match result {
            MnpResult::Ref(q) => format!("Ref({})", q),
            MnpResult::Alt(q) => format!("Alt({})", q),
            MnpResult::LowQuality => "LowQuality".to_string(),
            MnpResult::ThirdAllele => "ThirdAllele".to_string(),
            MnpResult::Structural => "Structural".to_string(),
        }
    }


    // ── SW vs PairHMM concordance tests ──
    //
    // These tests send identical synthetic reads through both backends
    // and assert that unambiguous cases produce the same REF/ALT decision.
    // Uses check_allele_with_qual (the real dispatch entry point).

    fn build_variant_with_context(
        pos: i64, ref_allele: &str, alt_allele: &str,
        ref_context: &str, ctx_start: i64,
    ) -> Variant {
        Variant {
            chrom: "1".to_string(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            variant_type: String::new(),
            ref_context: Some(ref_context.to_string()),
            ref_context_start: ctx_start,
            repeat_span: 0,
        }
    }

    /// Run check_allele_with_qual through both SW and HMM backends for concordance testing.
    fn run_both_backends(
        record: &Record, variant: &Variant, min_baseq: u8,
    ) -> ((bool, bool, u8), (bool, bool, u8)) {
        // We need separate aligner instances because check_allele_with_qual takes &mut
        let scoring_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let mut alt_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );
        let mut ref_a1 = Aligner::with_capacity_and_scoring(
            200, 200, bio::alignment::pairwise::Scoring::new(-5, -1, &scoring_fn)
                .xclip(bio::alignment::pairwise::MIN_SCORE)
                .yclip(0),
        );

        let sw_result = check_allele_with_qual(
            record, variant, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::SmithWaterman,
        );
        let hmm_result = check_allele_with_qual(
            record, variant, min_baseq,
            &mut alt_a1, &mut ref_a1,
            &AlignmentBackend::pairhmm_default(),
        );

        (sw_result, hmm_result)
    }

    #[test]
    fn test_concordance_snp_ref() {
        // Read matches REF haplotype at high quality — both backends should say REF
        let variant = build_variant_with_context(5, "C", "T", "GGGGGCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(11)]);
        let record = build_record(b"GGGGGCGGGGG", &[35_u8; 11], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert_eq!((sw.0, sw.1), (true, false), "SW should classify as REF");
        assert_eq!((hmm.0, hmm.1), (true, false), "PairHMM should classify as REF");
    }

    #[test]
    fn test_concordance_snp_alt() {
        // Read matches ALT haplotype at high quality — both backends should say ALT
        let variant = build_variant_with_context(5, "C", "T", "GGGGGCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(11)]);
        let record = build_record(b"GGGGGTGGGGG", &[35_u8; 11], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert_eq!((sw.0, sw.1), (false, true), "SW should classify as ALT");
        assert_eq!((hmm.0, hmm.1), (false, true), "PairHMM should classify as ALT");
    }

    #[test]
    fn test_concordance_insertion_alt() {
        // Read carries a 3bp insertion (A→ACCC) — both backends should agree ALT
        let variant = build_variant_with_context(4, "A", "ACCC", "GGGGAGGGGG", 0);
        let cigar = CigarString(vec![
            rust_htslib::bam::record::Cigar::Match(5),
            rust_htslib::bam::record::Cigar::Ins(3),
            rust_htslib::bam::record::Cigar::Match(5),
        ]);
        let record = build_record(b"GGGGACCCGGGGG", &[35_u8; 13], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert!(sw.1, "SW should classify insertion as ALT");
        assert!(hmm.1, "PairHMM should classify insertion as ALT");
    }

    #[test]
    fn test_concordance_deletion_alt() {
        // Read carries a 3bp deletion (ACCC→A) — both backends should agree ALT
        let variant = build_variant_with_context(4, "ACCC", "A", "GGGGACCCGGGGG", 0);
        let cigar = CigarString(vec![
            rust_htslib::bam::record::Cigar::Match(5),
            rust_htslib::bam::record::Cigar::Del(3),
            rust_htslib::bam::record::Cigar::Match(5),
        ]);
        let record = build_record(b"GGGGAGGGGG", &[35_u8; 10], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        assert!(sw.1, "SW should classify deletion as ALT");
        assert!(hmm.1, "PairHMM should classify deletion as ALT");
    }

    #[test]
    fn test_concordance_complex_delins() {
        // Complex variant: TC→GA (2bp substitution, same length) at pos 5
        // This is simpler than a DelIns — both haplotypes are the same length,
        // so both backends can classify reliably.
        // REF context: GGGGG TC GGGGG (13bp, variant at offset 5)
        // ALT read:    GGGGG GA GGGGG (matches ALT)
        let variant = build_variant_with_context(5, "TC", "GA", "GGGGGTCGGGGG", 0);
        let cigar = CigarString(vec![rust_htslib::bam::record::Cigar::Match(12)]);
        let record = build_record(b"GGGGGGAGGGGG", &[35_u8; 12], &cigar, 0);

        let (sw, hmm) = run_both_backends(&record, &variant, 20);
        // Both backends should agree on ALT for this unambiguous 2bp substitution
        assert!(sw.1, "SW should classify 2bp sub as ALT, got is_ref={} is_alt={}", sw.0, sw.1);
        assert!(hmm.1, "PairHMM should classify 2bp sub as ALT, got is_ref={} is_alt={}", hmm.0, hmm.1);
    }
}
