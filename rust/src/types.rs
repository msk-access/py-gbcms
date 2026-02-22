use pyo3::prelude::*;

#[pyclass]
#[derive(Debug, Clone)]
pub struct Variant {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub pos: i64, // 0-based
    #[pyo3(get, set)]
    pub ref_allele: String,
    #[pyo3(get, set)]
    pub alt_allele: String,
    #[pyo3(get, set)]
    pub variant_type: String, // "SNP", "INSERTION", "DELETION", "COMPLEX"
    /// Reference sequence around the variant for windowed indel detection.
    /// Covers [ref_context_start, ref_context_start + len) in genomic coords.
    /// Used by Safeguard 3 to verify shifted indels are biologically valid.
    #[pyo3(get, set)]
    pub ref_context: Option<String>,
    /// Genomic start position (0-based) of the ref_context string.
    #[pyo3(get, set)]
    pub ref_context_start: i64,
    /// Span of the tandem repeat region surrounding the variant (0 if not in a repeat).
    /// Used to dynamically tune Smith-Waterman gap penalties: repeat_span >= 10
    /// triggers gap_extend = 0 to absorb polymerase slippage noise.
    #[pyo3(get, set)]
    pub repeat_span: usize,
}

#[pymethods]
impl Variant {
    #[new]
    #[pyo3(signature = (chrom, pos, ref_allele, alt_allele, variant_type, ref_context=None, ref_context_start=0, repeat_span=0))]
    pub fn new(
        chrom: String,
        pos: i64,
        ref_allele: String,
        alt_allele: String,
        variant_type: String,
        ref_context: Option<String>,
        ref_context_start: i64,
        repeat_span: usize,
    ) -> Self {
        Variant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            variant_type,
            ref_context,
            ref_context_start,
            repeat_span,
        }
    }
}

#[pyclass]
#[derive(Debug, Clone, Default)]
pub struct BaseCounts {
    // Basic counts
    #[pyo3(get)]
    pub dp: u32,
    #[pyo3(get)]
    pub rd: u32,
    #[pyo3(get)]
    pub ad: u32,

    // Strand-specific counts
    #[pyo3(get)]
    pub dp_fwd: u32,
    #[pyo3(get)]
    pub rd_fwd: u32,
    #[pyo3(get)]
    pub ad_fwd: u32,
    #[pyo3(get)]
    pub dp_rev: u32,
    #[pyo3(get)]
    pub rd_rev: u32,
    #[pyo3(get)]
    pub ad_rev: u32,

    // Fragment counts (Majority Rule)
    #[pyo3(get)]
    pub dpf: u32,
    #[pyo3(get)]
    pub rdf: u32,
    #[pyo3(get)]
    pub adf: u32,

    // Fragment strand counts
    #[pyo3(get)]
    pub rdf_fwd: u32,
    #[pyo3(get)]
    pub rdf_rev: u32,
    #[pyo3(get)]
    pub adf_fwd: u32,
    #[pyo3(get)]
    pub adf_rev: u32,

    // Stats
    #[pyo3(get)]
    pub sb_pval: f64,
    #[pyo3(get)]
    pub sb_or: f64,
    #[pyo3(get)]
    pub fsb_pval: f64,
    #[pyo3(get)]
    pub fsb_or: f64,

    // Homopolymer decomposition
    /// True if the decomposed (corrected) allele was used because it
    /// produced a higher alt_count than the original allele.
    #[pyo3(get)]
    pub used_decomposed: bool,
}
