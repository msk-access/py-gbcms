# Type stubs for the Rust extension module (gbcms._rs)
# This file tells mypy about the types in the native extension

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
    chrom: str
    pos: int
    ref: str
    alt: str
    dp: int
    rd: int
    ad: int
    rd_fwd: int
    rd_rev: int
    ad_fwd: int
    ad_rev: int
    dp_fragment: int
    rd_fragment: int
    ad_fragment: int
    sb_pvalue: float

def count_bam(
    bam_path: str,
    variants: list[Variant],
    min_mapq: int = 20,
    min_baseq: int = 0,
    filter_duplicates: bool = True,
    filter_secondary: bool = False,
    filter_supplementary: bool = False,
    filter_qc_failed: bool = False,
    filter_improper_pair: bool = False,
    filter_indel: bool = False,
    threads: int = 1,
) -> list[BaseCounts]: ...
