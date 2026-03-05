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
    #[allow(clippy::too_many_arguments)]
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

    // ── mFSD: Fragment Size Distribution counts ───────────────────────────────
    // Number of fragments in each Krewlyzer size class (valid cfDNA range only;
    // 50–1000 bp with a non-zero TLEN). Gate downstream stats on mfsd_ks_valid.
    /// Fragments classified as REF and with a valid insert size.
    #[pyo3(get)]
    pub mfsd_ref_count: u32,
    /// Fragments classified as ALT and with a valid insert size.
    #[pyo3(get)]
    pub mfsd_alt_count: u32,
    /// Fragments that are neither REF nor ALT, nor an N base (third allele).
    #[pyo3(get)]
    pub mfsd_nonref_count: u32,
    /// Fragments where the base at the variant position was called 'N'.
    #[pyo3(get)]
    pub mfsd_n_count: u32,

    // ── mFSD: Mean fragment sizes ─────────────────────────────────────────────
    /// Mean fragment size (bp) for REF-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_ref_mean: f64,
    /// Mean fragment size (bp) for ALT-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_alt_mean: f64,
    /// Mean fragment size (bp) for NonREF-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_nonref_mean: f64,
    /// Mean fragment size (bp) for N-classified fragments. 0.0 when empty.
    #[pyo3(get)]
    pub mfsd_n_mean: f64,

    // ── mFSD: Log-Likelihood Ratios ───────────────────────────────────────────
    // LLR = Σ log(P_tumor(size) / P_healthy(size)) over all fragments in class.
    // Positive = tumor-like (short fragments); negative = healthy-like (long).
    /// LLR for ALT-classified fragments.
    #[pyo3(get)]
    pub mfsd_alt_llr: f64,
    /// LLR for REF-classified fragments.
    #[pyo3(get)]
    pub mfsd_ref_llr: f64,

    // ── mFSD: Pairwise KS comparisons (6 pairs × 3 values = 18 fields) ───────
    // Each triad: (delta = mean_A - mean_B, KS D-statistic, KS p-value).
    // NaN when either class has fewer than mfsd::MIN_FOR_KS (5) fragments.
    // Check mfsd_ks_valid (Python-derived) before interpreting these values.

    /// ALT vs REF: mean(ALT) − mean(REF)
    #[pyo3(get)]
    pub mfsd_delta_alt_ref: f64,
    /// ALT vs REF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_ref: f64,
    /// ALT vs REF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_ref: f64,

    /// ALT vs NonREF: mean(ALT) − mean(NonREF)
    #[pyo3(get)]
    pub mfsd_delta_alt_nonref: f64,
    /// ALT vs NonREF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_nonref: f64,
    /// ALT vs NonREF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_nonref: f64,

    /// REF vs NonREF: mean(REF) − mean(NonREF)
    #[pyo3(get)]
    pub mfsd_delta_ref_nonref: f64,
    /// REF vs NonREF: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_ref_nonref: f64,
    /// REF vs NonREF: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_ref_nonref: f64,

    /// ALT vs N: mean(ALT) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_alt_n: f64,
    /// ALT vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_alt_n: f64,
    /// ALT vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_alt_n: f64,

    /// REF vs N: mean(REF) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_ref_n: f64,
    /// REF vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_ref_n: f64,
    /// REF vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_ref_n: f64,

    /// NonREF vs N: mean(NonREF) − mean(N)
    #[pyo3(get)]
    pub mfsd_delta_nonref_n: f64,
    /// NonREF vs N: KS D-statistic
    #[pyo3(get)]
    pub mfsd_ks_nonref_n: f64,
    /// NonREF vs N: KS p-value
    #[pyo3(get)]
    pub mfsd_pval_nonref_n: f64,

    // ── mFSD: Raw size arrays (for --mfsd-parquet export) ────────────────────
    // Populated in all runs but only copied to disk when --mfsd-parquet is set.
    // NOT exported via PyO3 — written directly to Parquet by write_fsd_parquet()
    // in lib.rs, avoiding an FFI round-trip and the pyarrow dependency.
    /// Raw REF fragment sizes (bp). Internal only; use write_fsd_parquet() to persist.
    pub ref_sizes: Vec<u32>,
    /// Raw ALT fragment sizes (bp). Internal only; use write_fsd_parquet() to persist.
    pub alt_sizes: Vec<u32>,
}

