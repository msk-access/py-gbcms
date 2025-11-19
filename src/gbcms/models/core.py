"""
Core data models for gbcms v2.
"""

from enum import Enum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator


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
            raise ValueError(f"End position ({self.end}) must be >= start position ({self.start})")
        return self


class Variant(BaseModel):
    """
    Normalized representation of a genomic variant.
    """
    chrom: str
    pos: int = Field(ge=0, description="0-based position of the variant")
    ref: str
    alt: str
    variant_type: VariantType
    
    # Original input metadata (optional)
    original_id: str | None = None
    
    @property
    def interval(self) -> GenomicInterval:
        """
        Get the genomic interval covered by this variant.
        
        For SNP: [pos, pos+1)
        For Deletion: [pos, pos+len(ref)) (ref includes anchor? depends on normalization)
        For Insertion: [pos, pos+1) (anchor base)
        """
        # Note: This logic depends on strict normalization rules which we will implement in the kernel.
        # For now, a simple approximation based on ref length.
        return GenomicInterval(chrom=self.chrom, start=self.pos, end=self.pos + len(self.ref))


class InputFormat(str, Enum):
    VCF = "vcf"
    MAF = "maf"


class OutputFormat(str, Enum):
    VCF = "vcf"
    MAF = "maf"


class GbcmsConfig(BaseModel):
    """
    Global configuration for gbcms execution.
    """
    # Input
    variant_file: Path
    bam_files: dict[str, Path]  # sample_name -> bam_path
    reference_fasta: Path
    
    # Output
    output_dir: Path
    output_format: OutputFormat = OutputFormat.VCF
    
    # Filters
    min_mapping_quality: int = Field(default=20, ge=0)
    min_base_quality: int = Field(default=0, ge=0)
    filter_duplicates: bool = True
    filter_secondary: bool = False
    filter_supplementary: bool = False
    
    # Performance
    threads: int = Field(default=1, ge=1)
    
    # Advanced
    fragment_counting: bool = False
    
    @field_validator("variant_file", "reference_fasta")
    @classmethod
    def validate_file_exists(cls, v: Path) -> Path:
        if not v.exists():
            raise ValueError(f"File not found: {v}")
        return v

    @field_validator("output_dir")
    @classmethod
    def validate_output_dir(cls, v: Path) -> Path:
        if not v.exists():
            # Try to create it? Or just fail?
            # Usually safer to fail or let the pipeline create it.
            # But for config validation, let's just ensure it's not a file.
            if v.is_file():
                raise ValueError(f"Output path must be a directory, not a file: {v}")
        return v

    @model_validator(mode="after")
    def validate_bams(self) -> "GbcmsConfig":
        for name, path in self.bam_files.items():
            if not path.exists():
                raise ValueError(f"BAM file for sample '{name}' not found: {path}")
        return self
