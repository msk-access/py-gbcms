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
    """Quality score thresholds for filtering reads and bases."""

    min_mapping_quality: int = Field(default=20, ge=0, description="Minimum mapping quality (MAPQ)")
    min_base_quality: int = Field(default=0, ge=0, description="Minimum base quality (BQ)")
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

    @field_validator("directory")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        """Ensure output path is not a file."""
        if v.exists() and v.is_file():
            raise ValueError(f"Output path must be a directory, not a file: {v}")
        return v


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
