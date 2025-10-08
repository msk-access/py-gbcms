"""Tests for configuration module."""

from pathlib import Path

import pytest

from gbcms.config import Config, CountType


def test_count_type_enum():
    """Test CountType enum values."""
    assert CountType.DP.value == 0
    assert CountType.RD.value == 1
    assert CountType.AD.value == 2
    assert CountType.DPP.value == 3
    assert CountType.RDP.value == 4
    assert CountType.ADP.value == 5
    assert CountType.DPF.value == 6
    assert CountType.RDF.value == 7
    assert CountType.ADF.value == 8
    assert len(CountType) == 9


def test_config_validation_missing_fasta(temp_dir: Path, sample_bam: Path, sample_vcf: Path):
    """Test config validation with missing FASTA file."""
    with pytest.raises(FileNotFoundError, match="Reference FASTA file not found"):
        Config(
            fasta_file=str(temp_dir / "nonexistent.fa"),
            bam_files={"sample1": str(sample_bam)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_vcf=True,
        )


def test_config_validation_missing_fai(temp_dir: Path, sample_bam: Path, sample_vcf: Path):
    """Test config validation with missing FASTA index."""
    # Create FASTA without index (use different name to avoid fixture conflicts)
    fasta_file = temp_dir / "test_reference.fa"
    with open(fasta_file, "w") as f:
        f.write(">chr1\nATCG\n")

    with pytest.raises(FileNotFoundError, match="Reference FASTA index not found"):
        Config(
            fasta_file=str(fasta_file),
            bam_files={"sample1": str(sample_bam)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_vcf=True,
        )


def test_config_validation_missing_bam(temp_dir: Path, sample_fasta: Path, sample_vcf: Path):
    """Test config validation with missing BAM file."""
    with pytest.raises(FileNotFoundError, match="BAM file not found"):
        Config(
            fasta_file=str(sample_fasta),
            bam_files={"sample1": str(temp_dir / "nonexistent.bam")},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_vcf=True,
        )


def test_config_validation_missing_bai(temp_dir: Path, sample_fasta: Path, sample_vcf: Path):
    """Test config validation with missing BAM index."""
    # Create BAM without index
    bam_file = temp_dir / "sample.bam"
    bam_file.touch()

    with pytest.raises(FileNotFoundError, match="BAM index not found"):
        Config(
            fasta_file=str(sample_fasta),
            bam_files={"sample1": str(bam_file)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_vcf=True,
        )


def test_config_validation_mutually_exclusive(
    temp_dir: Path, sample_fasta: Path, sample_bam: Path, sample_vcf: Path
):
    """Test config validation with mutually exclusive options."""
    with pytest.raises(ValueError, match="mutually exclusive"):
        Config(
            fasta_file=str(sample_fasta),
            bam_files={"sample1": str(sample_bam)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_maf=True,
            input_is_vcf=True,
        )


def test_config_validation_missing_input_format(
    temp_dir: Path, sample_fasta: Path, sample_bam: Path, sample_vcf: Path
):
    """Test config validation without input format specified."""
    with pytest.raises(ValueError, match="Either --maf or --vcf must be specified"):
        Config(
            fasta_file=str(sample_fasta),
            bam_files={"sample1": str(sample_bam)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
        )


def test_config_validation_invalid_threads(
    temp_dir: Path, sample_fasta: Path, sample_bam: Path, sample_vcf: Path
):
    """Test config validation with invalid thread count."""
    with pytest.raises(ValueError, match="Number of threads must be at least 1"):
        Config(
            fasta_file=str(sample_fasta),
            bam_files={"sample1": str(sample_bam)},
            variant_files=[str(sample_vcf)],
            output_file=str(temp_dir / "output.txt"),
            input_is_vcf=True,
            num_threads=0,
        )


def test_config_valid(temp_dir: Path, sample_fasta: Path, sample_bam: Path, sample_vcf: Path):
    """Test valid configuration."""
    config = Config(
        fasta_file=str(sample_fasta),
        bam_files={"sample1": str(sample_bam)},
        variant_files=[str(sample_vcf)],
        output_file=str(temp_dir / "output.txt"),
        input_is_vcf=True,
        num_threads=4,
        mapping_quality_threshold=30,
        base_quality_threshold=20,
    )

    assert config.fasta_file == str(sample_fasta)
    assert config.num_threads == 4
    assert config.mapping_quality_threshold == 30
    assert config.base_quality_threshold == 20
    assert config.input_is_vcf is True
    assert config.input_is_maf is False
