//! Shared utility functions for the counting engine.
//!
//! Contains helper functions used by multiple counting submodules:
//! position lookup in BAM records, quality computation, masked
//! sequence comparison functions, and haplotype construction.

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use log::trace;

use crate::types::Variant;


/// Which classification phase resolved a read's allele assignment.
///
/// Used by `ClassifyResult` to track where in the multi-phase pipeline
/// each read was classified. Phases are ordered by computational cost
/// (cheapest first).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassifyPhase {
    /// Phase 0: Direct CIGAR pattern match (SNPs, structural indels)
    Structural,
    /// Phase 1: CIGAR-based sequence reconstruction + exact comparison
    CigarRecon,
    /// Phase 2: BQ-masked sequence comparison (tolerates low-quality bases)
    MaskedCompare,
    /// Phase 2.5: Levenshtein edit distance tiebreaker
    Levenshtein,
    /// Phase 3: Haplotype alignment (Smith-Waterman or PairHMM)
    Alignment,
}

/// Result of per-read allele classification.
///
/// Replaces the raw `(bool, bool, u8)` tuple to carry phase provenance,
/// enabling per-variant phase usage statistics.
#[derive(Debug, Clone, Copy)]
pub struct ClassifyResult {
    pub is_ref: bool,
    pub is_alt: bool,
    pub qual: u8,
    pub phase: ClassifyPhase,
}

impl ClassifyResult {
    /// Create a new ClassifyResult.
    #[inline]
    pub fn new(is_ref: bool, is_alt: bool, qual: u8, phase: ClassifyPhase) -> Self {
        Self { is_ref, is_alt, qual, phase }
    }

    /// Neither REF nor ALT — read didn't classify (e.g., no coverage, low quality).
    #[inline]
    pub fn neither(phase: ClassifyPhase) -> Self {
        Self { is_ref: false, is_alt: false, qual: 0, phase }
    }

    /// Shorthand for REF classification.
    #[inline]
    pub fn is_ref(qual: u8, phase: ClassifyPhase) -> Self {
        Self { is_ref: true, is_alt: false, qual, phase }
    }

    /// Shorthand for ALT classification.
    #[inline]
    pub fn is_alt(qual: u8, phase: ClassifyPhase) -> Self {
        Self { is_ref: false, is_alt: true, qual, phase }
    }

}

/// Build REF and ALT haplotypes from a variant's reference context.
///
/// Shared by both Smith-Waterman (`classify_by_alignment`) and PairHMM
/// (`classify_by_pairhmm`) backends. The haplotypes are constructed by
/// splicing the variant's alleles into the reference context:
///
///   ref_hap = left_ctx + REF_allele + right_ctx  (should equal ref_context)
///   alt_hap = left_ctx + ALT_allele + right_ctx
///
/// Returns `None` if the variant has no `ref_context` or the offset is invalid.
pub fn build_haplotypes(variant: &Variant) -> Option<(Vec<u8>, Vec<u8>)> {
    let ref_context = variant.ref_context.as_ref()?.as_bytes();
    let offset = (variant.pos - variant.ref_context_start) as usize;
    let ref_len = variant.ref_allele.len();

    if offset + ref_len > ref_context.len() {
        trace!(
            "build_haplotypes: offset {} + ref_len {} exceeds context len {}",
            offset, ref_len, ref_context.len()
        );
        return None;
    }

    let left_ctx = &ref_context[..offset];
    let right_ctx = &ref_context[offset + ref_len..];

    let ref_hap: Vec<u8> = left_ctx
        .iter()
        .chain(variant.ref_allele.as_bytes())
        .chain(right_ctx.iter())
        .copied()
        .collect();

    let alt_hap: Vec<u8> = left_ctx
        .iter()
        .chain(variant.alt_allele.as_bytes())
        .chain(right_ctx.iter())
        .copied()
        .collect();

    Some((ref_hap, alt_hap))
}

/// Find the read index corresponding to a genomic position.
pub fn find_read_pos(record: &Record, target_pos: i64) -> Option<usize> {
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

/// Compute the median quality of bases that pass the minimum threshold.
///
/// Follows the GATK `BaseQuality` annotation standard (median rather than min)
/// to prevent a single low-quality outlier from penalizing an entire read's
/// contribution to fragment consensus.
///
/// Returns 0 if no qualifying bases.
#[inline]
pub fn median_qual(quals: &[u8], min_baseq: u8) -> u8 {
    let mut filtered: Vec<u8> = quals.iter()
        .copied()
        .filter(|&q| q >= min_baseq)
        .collect();
    if filtered.is_empty() { return 0; }
    filtered.sort_unstable();
    filtered[filtered.len() / 2]
}

/// Masked comparison against two alleles simultaneously.
///
/// Masks out bases with quality below `min_baseq` and counts mismatches against
/// both `allele_a` and `allele_b` on the remaining reliable bases only.
///
/// Returns `(mismatches_a, mismatches_b, reliable_count)`.
pub fn masked_dual_compare(
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
pub fn masked_single_compare(
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
