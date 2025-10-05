"""Configuration classes and enums for GetBaseCounts."""

import os
from dataclasses import dataclass
from enum import IntEnum


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


@dataclass
class Config:
    """Configuration for base counting."""

    fasta_file: str
    bam_files: dict[str, str]  # sample_name -> bam_path
    variant_files: list[str]
    output_file: str

    # Optional parameters
    mapping_quality_threshold: int = 20
    base_quality_threshold: int = 0
    filter_duplicate: bool = True
    filter_improper_pair: bool = False
    filter_qc_failed: bool = False
    filter_indel: bool = False
    filter_non_primary: bool = False
    output_positive_count: bool = True
    output_negative_count: bool = False
    output_fragment_count: bool = False
    fragment_fractional_weight: bool = False
    max_block_size: int = 10000
    max_block_dist: int = 100000
    num_threads: int = 1
    backend: str = "joblib"  # Parallelization backend
    input_is_maf: bool = False
    input_is_vcf: bool = False
    output_maf: bool = False
    generic_counting: bool = False
    max_warning_per_type: int = 3

    def __post_init__(self) -> None:
        """Validate configuration."""
        if not os.path.exists(self.fasta_file):
            raise FileNotFoundError(f"Reference FASTA file not found: {self.fasta_file}")

        fai_file = f"{self.fasta_file}.fai"
        if not os.path.exists(fai_file):
            raise FileNotFoundError(
                f"Reference FASTA index not found: {fai_file}. "
                f"Please index with: samtools faidx {self.fasta_file}"
            )

        for sample, bam_path in self.bam_files.items():
            if not os.path.exists(bam_path):
                raise FileNotFoundError(f"BAM file not found for sample {sample}: {bam_path}")

            # Check for BAM index
            bai_file1 = bam_path.replace(".bam", ".bai")
            bai_file2 = f"{bam_path}.bai"
            if not os.path.exists(bai_file1) and not os.path.exists(bai_file2):
                raise FileNotFoundError(
                    f"BAM index not found for {bam_path}. "
                    f"Please index with: samtools index {bam_path}"
                )

        for variant_file in self.variant_files:
            if not os.path.exists(variant_file):
                raise FileNotFoundError(f"Variant file not found: {variant_file}")

        if self.input_is_maf and self.input_is_vcf:
            raise ValueError("--maf and --vcf are mutually exclusive")

        if not self.input_is_maf and not self.input_is_vcf:
            raise ValueError("Either --maf or --vcf must be specified")

        if self.input_is_vcf and self.output_maf:
            raise ValueError("--omaf can only be used with --maf input")

        if self.num_threads < 1:
            raise ValueError("Number of threads must be at least 1")

        if self.max_block_size < 1:
            raise ValueError("max_block_size must be at least 1")

        if self.max_block_dist < 1:
            raise ValueError("max_block_dist must be at least 1")
