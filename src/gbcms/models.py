"""Pydantic models for type-safe configuration and data structures."""

from enum import IntEnum
from pathlib import Path

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator


class CountType(IntEnum):
    """Enumeration for different count types."""

    DP = 0  # Total depth
    RD = 1  # Reference depth
    AD = 2  # Alternate depth
    DPP = 3  # Positive strand depth
    RDP = 4  # Positive strand reference depth
    ADP = 5  # Positive strand alternate depth
    DPF = 6  # Fragment depth
    RDF = 7  # Fragment reference depth
    ADF = 8  # Fragment alternate depth


class BamFileConfig(BaseModel):
    """Configuration for a single BAM file."""

    sample_name: str = Field(..., description="Sample name")
    bam_path: Path = Field(..., description="Path to BAM file")
    bai_path: Path | None = Field(None, description="Path to BAM index")

    @field_validator("bam_path")
    @classmethod
    def validate_bam_exists(cls, v: Path) -> Path:
        """Validate BAM file exists."""
        if not v.exists():
            raise ValueError(f"BAM file not found: {v}")
        return v

    @model_validator(mode="after")
    def validate_bai(self) -> "BamFileConfig":
        """Validate BAM index exists."""
        if self.bai_path is None:
            # Try to find index
            bai_path1 = Path(str(self.bam_path).replace(".bam", ".bai"))
            bai_path2 = Path(f"{self.bam_path}.bai")

            if bai_path1.exists():
                self.bai_path = bai_path1
            elif bai_path2.exists():
                self.bai_path = bai_path2
            else:
                raise ValueError(f"BAM index not found for: {self.bam_path}")

        return self

    model_config = {"arbitrary_types_allowed": True}


class VariantFileConfig(BaseModel):
    """Configuration for variant files."""

    file_path: Path = Field(..., description="Path to variant file")
    file_format: str = Field(..., description="File format (vcf or maf)")

    @field_validator("file_path")
    @classmethod
    def validate_file_exists(cls, v: Path) -> Path:
        """Validate variant file exists."""
        if not v.exists():
            raise ValueError(f"Variant file not found: {v}")
        return v

    @field_validator("file_format")
    @classmethod
    def validate_format(cls, v: str) -> str:
        """Validate file format."""
        if v.lower() not in ["vcf", "maf"]:
            raise ValueError(f"Invalid format: {v}. Must be 'vcf' or 'maf'")
        return v.lower()

    model_config = {"arbitrary_types_allowed": True}


class QualityFilters(BaseModel):
    """Quality filtering parameters."""

    mapping_quality_threshold: int = Field(20, ge=0, description="Mapping quality threshold")
    base_quality_threshold: int = Field(0, ge=0, description="Base quality threshold")
    filter_duplicate: bool = Field(True, description="Filter duplicate reads")
    filter_improper_pair: bool = Field(False, description="Filter improper pairs")
    filter_qc_failed: bool = Field(False, description="Filter QC failed reads")
    filter_indel: bool = Field(False, description="Filter reads with indels")
    filter_non_primary: bool = Field(False, description="Filter non-primary alignments")


class OutputOptions(BaseModel):
    """Output configuration options."""

    output_file: Path = Field(..., description="Output file path")
    output_maf: bool = Field(False, description="Output in MAF format")
    output_positive_count: bool = Field(True, description="Output positive strand counts")
    output_negative_count: bool = Field(False, description="Output negative strand counts")
    output_fragment_count: bool = Field(False, description="Output fragment counts")
    fragment_fractional_weight: bool = Field(
        False, description="Use fractional weights for fragments"
    )

    model_config = {"arbitrary_types_allowed": True}


class PerformanceConfig(BaseModel):
    """Performance and parallelization configuration."""

    num_threads: int = Field(1, ge=1, description="Number of threads")
    max_block_size: int = Field(10000, ge=1, description="Maximum variants per block")
    max_block_dist: int = Field(100000, ge=1, description="Maximum block distance in bp")
    use_numba: bool = Field(True, description="Use Numba JIT compilation")

    @field_validator("backend")
    @classmethod
    def validate_backend(cls, v: str) -> str:
        """Validate backend choice."""
        valid_backends = ["joblib", "loky", "threading", "multiprocessing"]
        if v.lower() not in valid_backends:
            raise ValueError(f"Invalid backend: {v}. Must be one of: {', '.join(valid_backends)}")
        return v.lower()


class GetBaseCountsConfig(BaseModel):
    """Complete configuration for GetBaseCounts with Pydantic validation."""

    # Input files
    fasta_file: Path = Field(..., description="Reference FASTA file")
    bam_files: list[BamFileConfig] = Field(..., description="BAM files to process")
    variant_files: list[VariantFileConfig] = Field(..., description="Variant files")

    # Options
    quality_filters: QualityFilters = Field(
        default_factory=QualityFilters, description="Quality filtering options"  # type: ignore[arg-type]
    )
    output_options: OutputOptions = Field(..., description="Output options")
    performance: PerformanceConfig = Field(
        default_factory=PerformanceConfig, description="Performance options"  # type: ignore[arg-type]
    )

    # Advanced
    generic_counting: bool = Field(False, description="Use generic counting algorithm")
    max_warning_per_type: int = Field(3, ge=0, description="Maximum warnings per type")

    @field_validator("fasta_file")
    @classmethod
    def validate_fasta_exists(cls, v: Path) -> Path:
        """Validate FASTA file exists."""
        if not v.exists():
            raise ValueError(f"FASTA file not found: {v}")

        fai_file = Path(f"{v}.fai")
        if not fai_file.exists():
            raise ValueError(f"FASTA index not found: {fai_file}")

        return v

    @model_validator(mode="after")
    def validate_variant_format_consistency(self) -> "GetBaseCountsConfig":
        """Validate variant file format consistency."""
        formats = {vf.file_format for vf in self.variant_files}
        if len(formats) > 1:
            raise ValueError("All variant files must be the same format (all VCF or all MAF)")

        # Check MAF output compatibility
        if self.output_options.output_maf and "maf" not in formats:
            raise ValueError("--omaf can only be used with MAF input")

        return self

    def get_sample_names(self) -> list[str]:
        """Get list of sample names in order."""
        return [bam.sample_name for bam in self.bam_files]

    def is_maf_input(self) -> bool:
        """Check if input is MAF format."""
        return self.variant_files[0].file_format == "maf"

    def is_vcf_input(self) -> bool:
        """Check if input is VCF format."""
        return self.variant_files[0].file_format == "vcf"

    model_config = {"arbitrary_types_allowed": True}


class VariantCounts(BaseModel):
    """Type-safe variant counts structure."""

    sample_name: str
    counts: np.ndarray = Field(..., description="Count array")

    @field_validator("counts")
    @classmethod
    def validate_counts_shape(cls, v: np.ndarray) -> np.ndarray:
        """Validate counts array shape."""
        if v.shape != (len(CountType),):
            raise ValueError(f"Counts array must have shape ({len(CountType)},)")
        return v

    def get_count(self, count_type: CountType) -> float:
        """Get count for specific type."""
        return float(self.counts[count_type])

    def set_count(self, count_type: CountType, value: float) -> None:
        """Set count for specific type."""
        self.counts[count_type] = value

    model_config = {"arbitrary_types_allowed": True}


class VariantModel(BaseModel):
    """Pydantic model for variant with type safety."""

    chrom: str = Field(..., description="Chromosome")
    pos: int = Field(..., ge=0, description="Position (0-indexed)")
    end_pos: int = Field(..., ge=0, description="End position")
    ref: str = Field(..., min_length=1, description="Reference allele")
    alt: str = Field(..., min_length=1, description="Alternate allele")

    # Variant type flags
    snp: bool = Field(False, description="Is SNP")
    dnp: bool = Field(False, description="Is DNP")
    dnp_len: int = Field(0, ge=0, description="DNP length")
    insertion: bool = Field(False, description="Is insertion")
    deletion: bool = Field(False, description="Is deletion")

    # Sample information
    tumor_sample: str = Field("", description="Tumor sample name")
    normal_sample: str = Field("", description="Normal sample name")

    # Annotation
    gene: str = Field("", description="Gene name")
    effect: str = Field("", description="Variant effect")

    # Original MAF coordinates
    maf_pos: int = Field(0, ge=0, description="Original MAF position")
    maf_end_pos: int = Field(0, ge=0, description="Original MAF end position")
    maf_ref: str = Field("", description="Original MAF reference")
    maf_alt: str = Field("", description="Original MAF alternate")
    caller: str = Field("", description="Variant caller")

    # Counts
    sample_counts: dict[str, VariantCounts] = Field(
        default_factory=dict, description="Counts per sample"
    )

    @model_validator(mode="after")
    def validate_positions(self) -> "VariantModel":
        """Validate position consistency."""
        if self.end_pos < self.pos:
            raise ValueError(f"End position {self.end_pos} < start position {self.pos}")
        return self

    @model_validator(mode="after")
    def validate_variant_type(self) -> "VariantModel":
        """Validate variant type flags are consistent."""
        type_count = sum([self.snp, self.dnp, self.insertion, self.deletion])
        if type_count == 0:
            # Auto-detect variant type
            if len(self.ref) == len(self.alt) == 1:
                self.snp = True
            elif len(self.ref) == len(self.alt) > 1:
                self.dnp = True
                self.dnp_len = len(self.ref)
            elif len(self.alt) > len(self.ref):
                self.insertion = True
            elif len(self.alt) < len(self.ref):
                self.deletion = True

        return self

    def get_variant_key(self) -> tuple[str, int, str, str]:
        """Get unique variant key."""
        return (self.chrom, self.pos, self.ref, self.alt)

    def initialize_counts(self, sample_names: list[str]) -> None:
        """Initialize counts for all samples."""
        for sample in sample_names:
            if sample not in self.sample_counts:
                self.sample_counts[sample] = VariantCounts(
                    sample_name=sample, counts=np.zeros(len(CountType), dtype=np.float32)
                )

    def get_count(self, sample: str, count_type: CountType) -> float:
        """Get count for specific sample and type."""
        if sample not in self.sample_counts:
            return 0.0
        return self.sample_counts[sample].get_count(count_type)

    model_config = {"arbitrary_types_allowed": True}
