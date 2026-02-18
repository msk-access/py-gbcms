# Type stubs for the Rust extension module (gbcms._rs)
# This file tells mypy about the types in the native extension

class Variant:
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    variant_type: str
    ref_context: str | None
    ref_context_start: int

    def __init__(
        self,
        chrom: str,
        pos: int,
        ref_allele: str,
        alt_allele: str,
        variant_type: str,
        ref_context: str | None = None,
        ref_context_start: int = 0,
    ) -> None: ...

class BaseCounts:
    chrom: str
    pos: int
    ref: str
    alt: str
    dp: int
    rd: int
    ad: int
    dp_fwd: int
    rd_fwd: int
    ad_fwd: int
    dp_rev: int
    rd_rev: int
    ad_rev: int
    dpf: int
    rdf: int
    adf: int
    rdf_fwd: int
    rdf_rev: int
    adf_fwd: int
    adf_rev: int
    sb_pval: float
    sb_or: float
    fsb_pval: float
    fsb_or: float
    used_decomposed: bool

class PreparedVariant:
    variant: Variant
    validation_status: str
    was_normalized: bool
    original_pos: int
    original_ref: str
    original_alt: str
    decomposed_variant: Variant | None

def count_bam(
    bam_path: str,
    variants: list[Variant],
    decomposed: list[Variant | None],
    min_mapq: int = 20,
    min_baseq: int = 20,
    filter_duplicates: bool = True,
    filter_secondary: bool = False,
    filter_supplementary: bool = False,
    filter_qc_failed: bool = False,
    filter_improper_pair: bool = False,
    filter_indel: bool = False,
    threads: int = 1,
    fragment_qual_threshold: int = 10,
) -> list[BaseCounts]: ...
def prepare_variants(
    variants: list[Variant],
    fasta_path: str,
    context_padding: int,
    is_maf: bool,
    threads: int = 1,
    adaptive_context: bool = True,
) -> list[PreparedVariant]: ...
