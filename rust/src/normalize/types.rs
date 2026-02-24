//! `PreparedVariant` — output type from variant preparation.

use pyo3::prelude::*;
use crate::types::Variant;

/// Result of variant preparation: normalized coords, ref_context, and validation info.
///
/// Every input variant produces exactly one `PreparedVariant`, even if validation
/// fails — this ensures the output always has the same row count as input.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PreparedVariant {
    /// Ready-to-count variant (normalized coords + ref_context populated).
    /// For invalid variants, contains best-effort coords with no ref_context.
    #[pyo3(get)]
    pub variant: Variant,

    /// Validation result: "PASS", "PASS_WARN_REF_CORRECTED",
    /// "PASS_WARN_HOMOPOLYMER_DECOMP", "REF_MISMATCH", or "FETCH_FAILED".
    #[pyo3(get, set)]
    pub validation_status: String,

    /// True if left-alignment changed the variant's coordinates.
    #[pyo3(get)]
    pub was_normalized: bool,

    /// Original 0-based position before any transformation.
    #[pyo3(get)]
    pub original_pos: i64,

    /// Original REF allele before any transformation.
    #[pyo3(get)]
    pub original_ref: String,

    /// Original ALT allele before any transformation.
    #[pyo3(get)]
    pub original_alt: String,

    /// Corrected variant for homopolymer decomposition dual-counting.
    /// When a complex variant spans a homopolymer and appears to be a
    /// miscollapsed D(n)+SNV event, this holds the corrected allele
    /// (e.g., CCCCCC→CCCCT instead of CCCCCC→T).
    /// `None` for normal variants where no decomposition is detected.
    #[pyo3(get)]
    pub decomposed_variant: Option<Variant>,

    /// Group ID for overlapping multi-allelic variants at the same locus.
    /// `None` for isolated variants, `Some(id)` when multiple variants share
    /// overlapping genomic footprints (same chrom, overlapping REF spans).
    #[pyo3(get)]
    pub multi_allelic_group: Option<u32>,
}
