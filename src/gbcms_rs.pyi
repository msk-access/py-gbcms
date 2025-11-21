class Variant:
    chrom: str
    pos: int
    ref_allele: str
    alt_allele: str
    variant_type: str
    def __init__(
        self,
        chrom: str,
        pos: int,
        ref_allele: str,
        alt_allele: str,
        variant_type: str,
    ) -> None: ...

class BaseCounts:
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

def count_bam(
    bam_path: str,
    variants: list[Variant],
    min_mapq: int,
    min_baseq: int,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: int,
) -> list[BaseCounts]: ...
