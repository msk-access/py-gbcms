"""
Core data models for gbcms v2.

This module defines the data models for variants, configuration, and nested
config groups (filters, quality thresholds, output settings).
"""

from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator, model_validator

__all__ = [
    "VariantType",
    "GenomicInterval",
    "Variant",
    "OutputFormat",
    "ReadFilters",
    "QualityThresholds",
    "OutputConfig",
    "GbcmsConfig",
]


class VariantType(str, Enum):
    """Type of genomic variant."""

    SNP = "SNP"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    COMPLEX = "COMPLEX"


class GenomicInterval(BaseModel):
    """
    Represents a 0-based, half-open genomic interval [start, end).

    This is the canonical internal representation for all coordinates.
    """

    chrom: str
    start: int = Field(ge=0, description="0-based start position (inclusive)")
    end: int = Field(ge=0, description="0-based end position (exclusive)")

    @model_validator(mode="after")
    def validate_interval(self) -> "GenomicInterval":
        if self.end < self.start:
            raise ValueError(
                f"End position ({self.end}) must be >= start position ({self.start})"
            )
        return self


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

    @property
    def interval(self) -> GenomicInterval:
        """Get the genomic interval covered by this variant."""
        return GenomicInterval(
            chrom=self.chrom, start=self.pos, end=self.pos + len(self.ref)
        )


class OutputFormat(str, Enum):
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
    supplementary: bool = Field(
        default=False, description="Filter supplementary alignments"
    )
    qc_failed: bool = Field(default=False, description="Filter reads failing QC")
    improper_pair: bool = Field(
        default=False, description="Filter improperly paired reads"
    )
    indel: bool = Field(default=False, description="Filter reads containing indels")


class QualityThresholds(BaseModel):
    """Quality score thresholds for filtering reads and bases."""

    min_mapping_quality: int = Field(
        default=20, ge=0, description="Minimum mapping quality (MAPQ)"
    )
    min_base_quality: int = Field(
        default=0, ge=0, description="Minimum base quality (BQ)"
    )


class OutputConfig(BaseModel):
    """Output configuration settings."""

    directory: Path = Field(description="Directory to write output files")
    format: OutputFormat = Field(
        default=OutputFormat.VCF, description="Output format (vcf or maf)"
    )
    suffix: str = Field(
        default="", description="Suffix to append to output filename"
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
    fragment_counting: bool = Field(
        default=False, description="Enable fragment-based counting"
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
