"""
Core data models for gbcms v2.

This module defines the data models for variants, configuration, and nested
config groups (filters, quality thresholds, output settings).
"""

import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from enum import Enum

    class StrEnum(str, Enum):
        """Backport of StrEnum for Python 3.10."""

        pass


from pydantic import BaseModel, Field, field_validator, model_validator

__all__ = [
    "VariantType",
    "Variant",
    "OutputFormat",
    "ReadFilters",
    "QualityThresholds",
    "AlignmentConfig",
    "OutputConfig",
    "GbcmsConfig",
]


class VariantType(StrEnum):
    """Type of genomic variant."""

    SNP = "SNP"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    COMPLEX = "COMPLEX"


class Variant(BaseModel):
    """Normalized representation of a genomic variant."""

    chrom: str
    pos: int = Field(ge=0, description="0-based position of the variant")
    ref: str
    alt: str
    variant_type: VariantType

    # Original input metadata (optional)
    original_id: str | None = None
    metadata: dict[str, str] = Field(
        default_factory=dict, description="Original input metadata/columns"
    )


class OutputFormat(StrEnum):
    """Supported output formats for gbcms."""

    VCF = "vcf"
    MAF = "maf"


# =============================================================================
# Nested Configuration Models
# =============================================================================


class ReadFilters(BaseModel):
    """
    Filters for read selection during BAM processing.

    These flags control which reads are excluded from counting.
    When True, reads with the corresponding flag are filtered out.
    """

    duplicates: bool = Field(default=True, description="Filter duplicate reads")
    secondary: bool = Field(default=False, description="Filter secondary alignments")
    supplementary: bool = Field(default=False, description="Filter supplementary alignments")
    qc_failed: bool = Field(default=False, description="Filter reads failing QC")
    improper_pair: bool = Field(default=False, description="Filter improperly paired reads")
    indel: bool = Field(default=False, description="Filter reads containing indels")


class QualityThresholds(BaseModel):
    """Quality score thresholds for filtering reads and bases.

    Defaults are chosen to match the CLI defaults so that programmatic callers
    and CLI users observe identical behaviour without explicit overrides.
    """

    min_mapping_quality: int = Field(
        default=20, ge=0, description="Minimum mapping quality (MAPQ). Default 20 matches CLI."
    )
    min_base_quality: int = Field(
        default=20,
        ge=0,
        description=(
            "Minimum base quality (BQ). Default 20 matches CLI `--min-baseq` default. "
            "Previously this was 0 (maximally permissive), which diverged from the CLI."
        ),
    )
    fragment_qual_threshold: int = Field(
        default=10,
        ge=0,
        le=93,
        description=(
            "Quality difference threshold for fragment consensus. "
            "When R1 and R2 disagree, the allele with higher base quality wins "
            "only if the difference exceeds this threshold. Ambiguous fragments "
            "(within threshold) are discarded to preserve VAF accuracy."
        ),
    )
    context_padding: int = Field(
        default=5,
        ge=1,
        le=50,
        description=(
            "Number of flanking reference bases fetched around indel/complex "
            "variants for haplotype construction and Smith-Waterman alignment. "
            "Larger values improve sensitivity for shifted indels in repeat "
            "regions at minimal computational cost."
        ),
    )
    adaptive_context: bool = Field(
        default=True,
        description=(
            "Dynamically increase context padding in tandem repeat regions. "
            "When enabled, the effective padding is max(context_padding, "
            "repeat_span/2 + 3), capped at 50bp."
        ),
    )


class AlignmentConfig(BaseModel):
    """Alignment backend configuration for Phase 3 read classification.

    Controls which algorithm is used for haplotype-level alignment when
    variant-type-specific CIGAR matching is inconclusive. Smith-Waterman (sw)
    is the default; PairHMM (hmm) is an alternative that uses per-base
    quality-aware emission probabilities.
    """

    backend: str = Field(
        default="sw",
        description="Alignment backend: 'sw' (Smith-Waterman) or 'hmm' (PairHMM).",
    )
    hmm_llr_threshold: float = Field(
        default=2.3,
        gt=0.0,
        description=(
            "PairHMM log-likelihood ratio threshold for confident calls. "
            "ln(10) ≈ 2.3 = 10:1 odds ratio (default)."
        ),
    )
    hmm_gap_open: float = Field(
        default=1e-4,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-open probability for non-repeat regions.",
    )
    hmm_gap_extend: float = Field(
        default=0.1,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-extend probability for non-repeat regions.",
    )
    hmm_gap_open_repeat: float = Field(
        default=1e-2,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-open probability for tandem repeat regions.",
    )
    hmm_gap_extend_repeat: float = Field(
        default=0.5,
        gt=0.0,
        lt=1.0,
        description="PairHMM gap-extend probability for tandem repeat regions.",
    )

    @field_validator("backend")
    @classmethod
    def validate_backend(cls, v: str) -> str:
        """Validate backend is a supported value."""
        v = v.lower().strip()
        if v not in ("sw", "hmm", "pairhmm"):
            raise ValueError(f"Invalid alignment backend '{v}'. Must be 'sw' or 'hmm'.")
        return v


class OutputConfig(BaseModel):
    """Output configuration settings."""

    directory: Path = Field(description="Directory to write output files")
    format: OutputFormat = Field(default=OutputFormat.VCF, description="Output format (vcf or maf)")
    suffix: str = Field(default="", description="Suffix to append to output filename")
    column_prefix: str = Field(
        default="",
        description=(
            "Prefix for gbcms count columns in MAF output. "
            "Default: no prefix (e.g., 'ref_count'). "
            "Use 't_' for legacy compatibility (e.g., 't_ref_count')."
        ),
    )
    preserve_barcode: bool = Field(
        default=False,
        description=(
            "When True, preserve the original Tumor_Sample_Barcode from "
            "input MAF instead of overriding with the BAM sample name. "
            "Only applies to MAF→MAF output; VCF→MAF always uses BAM name."
        ),
    )
    mfsd: bool = Field(
        default=False,
        description=(
            "Enable Mutant Fragment Size Distribution (mFSD) analysis. "
            "Adds 31 mFSD columns (KS test, LLR, mean sizes, pairwise comparisons, "
            "derived metrics) to MAF output and 7 MFSD INFO fields to VCF. "
            "Required when mfsd_parquet=True."
        ),
    )
    mfsd_parquet: bool = Field(
        default=False,
        description=(
            "Write a companion <sample>.fsd.parquet file with per-variant raw "
            "fragment size arrays (REF and ALT insert sizes). "
            "Enables downstream mFSD visualizations (e.g. density plots). "
            "Requires mfsd=True."
        ),
    )

    @field_validator("directory")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Ensure output path is not a file."""
        if v.exists() and v.is_file():
            raise ValueError(f"Output path must be a directory, not a file: {v}")
        return v

    @model_validator(mode="after")
    def validate_mfsd_parquet(self) -> "OutputConfig":
        """Enforce that mfsd_parquet requires mfsd to be enabled.

        This is validated at model construction so both CLI and programmatic
        callers get the same fail-fast behaviour.
        """
        if self.mfsd_parquet and not self.mfsd:
            raise ValueError(
                "mfsd_parquet=True requires mfsd=True. "
                "Enable mFSD analysis before requesting Parquet export."
            )
        return self


class GbcmsConfig(BaseModel):
    """
    Global configuration for gbcms execution.

    Groups related settings into nested models for cleaner organization.
    """

    # Input files
    variant_file: Path
    bam_files: dict[str, Path]  # sample_name -> bam_path
    reference_fasta: Path

    # Nested configuration groups
    output: OutputConfig
    filters: ReadFilters = Field(default_factory=ReadFilters)
    quality: QualityThresholds = Field(default_factory=QualityThresholds)

    # Performance
    threads: int = Field(default=1, ge=1, description="Number of threads")

    # Alignment backend
    alignment: AlignmentConfig = Field(default_factory=AlignmentConfig)

    # Advanced
    show_normalization: bool = Field(
        default=False,
        description="Add normalization columns showing left-aligned coordinates to output.",
    )

    @field_validator("variant_file", "reference_fasta")
    @classmethod
    def validate_file_exists(cls, v: Path) -> Path:
        """Validate that input files exist."""
        if not v.exists():
            raise ValueError(f"File not found: {v}")
        return v

    @model_validator(mode="after")
    def validate_bams(self) -> "GbcmsConfig":
        """Validate that all BAM files exist."""
        for name, path in self.bam_files.items():
            if not path.exists():
                raise ValueError(f"BAM file for sample '{name}' not found: {path}")
        return self
