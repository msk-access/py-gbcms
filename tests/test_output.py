"""Tests for output module."""

from pathlib import Path

import pytest

from gbcms.config import Config, CountType
from gbcms.output import SampleAgnosticMAFWriter, VCFWriter
from gbcms.variant import VariantEntry


@pytest.fixture
def config_vcf(temp_dir, sample_fasta, sample_bam, sample_vcf):
    """Create a VCF config."""
    return Config(
        fasta_file=str(sample_fasta),
        bam_files={"sample1": str(sample_bam), "sample2": str(sample_bam)},
        variant_files=[str(sample_vcf)],
        output_file=str(temp_dir / "output.txt"),
        input_is_vcf=True,
    )


@pytest.fixture
def config_maf(temp_dir, sample_fasta, sample_bam, sample_maf):
    """Create a MAF config."""
    return Config(
        fasta_file=str(sample_fasta),
        bam_files={"Tumor1": str(sample_bam), "Normal1": str(sample_bam)},
        variant_files=[str(sample_maf)],
        output_file=str(temp_dir / "output.maf"),
        input_is_maf=True,
    )


def test_vcf_writer(config_vcf, temp_dir):
    """Test VCFWriter formatting."""
    writer = VCFWriter(config_vcf, ["sample1", "sample2"])

    # Create test variants
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )
    variant.initialize_counts(["sample1", "sample2"])
    variant.base_count["sample1"][CountType.DP] = 10
    variant.base_count["sample1"][CountType.RD] = 6
    variant.base_count["sample1"][CountType.AD] = 4

    writer.write_variants([variant])

    # Check output file exists
    assert Path(config_vcf.output_file).exists()

    # Read and verify content (proper VCF format)
    with open(config_vcf.output_file) as f:
        lines = f.readlines()
        # Should have header lines + 1 variant line
        assert len(lines) >= 3  # Multiple header lines + 1 variant
        # Check for proper VCF headers
        content = "\n".join(lines)
        assert "##fileformat=VCFv4.2" in content
        assert "#CHROM" in content
        assert "chr1" in content


def test_sample_agnostic_maf_writer(temp_dir):
    """Test SampleAgnosticMAFWriter formatting."""
    # Create a simple config for testing
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={
            "Tumor1": "tests/testdata/sample1_integration_test.bam",
            "Normal1": "tests/testdata/sample2_integration_test.bam",
        },
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file=str(temp_dir / "output.maf"),
        input_is_maf=True,
    )

    writer = SampleAgnosticMAFWriter(config, ["Tumor1", "Normal1"])

    # Create test variant
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
        gene="GENE1",
        tumor_sample="Tumor1",
        normal_sample="Normal1",
        effect="Missense_Mutation",
    )
    variant.initialize_counts(["Tumor1", "Normal1"])
    variant.base_count["Tumor1"][CountType.DP] = 10
    variant.base_count["Tumor1"][CountType.RD] = 6
    variant.base_count["Tumor1"][CountType.AD] = 4

    # Write variants
    writer.write_variants([variant])

    # Check output file exists and has content
    assert Path(config.output_file).exists()

    with open(config.output_file) as f:
        lines = f.readlines()
        assert len(lines) >= 2  # Header + at least one variant line

        # Check for MAF header
        header_line = lines[0]
        assert "Hugo_Symbol" in header_line
        assert "Chromosome" in header_line

        # Check for variant data
        variant_line = lines[1]
        assert "GENE1" in variant_line
        assert "chr1" in variant_line
        assert "Missense_Mutation" in variant_line


def test_output_formatter_fillout(temp_dir):
    """Test fillout output formatting."""
    # Create a simple config for testing
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"Tumor1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file=str(temp_dir / "output.maf"),
        input_is_maf=True,
    )
    writer = SampleAgnosticMAFWriter(config, ["Tumor1"])

    # Create test variant
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
        gene="GENE1",
        tumor_sample="Tumor1",
        normal_sample="Normal1",
        effect="Missense_Mutation",
    )
    variant.initialize_counts(["Tumor1", "Normal1"])

    writer.write_variants([variant])

    # Check output file exists
    assert Path(config.output_file).exists()


def test_output_formatter_with_strand_counts(config_vcf, temp_dir):
    """Test output with strand counts."""
    config_vcf.output_positive_count = True
    writer = VCFWriter(config_vcf, ["sample1"])

    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )
    variant.initialize_counts(["sample1"])
    variant.base_count["sample1"][CountType.DP] = 10
    variant.base_count["sample1"][CountType.RD] = 6
    variant.base_count["sample1"][CountType.AD] = 4

    writer.write_variants([variant])

    # Check output includes strand counts in FORMAT field
    with open(config_vcf.output_file) as f:
        content = f.read()
        assert "DP" in content  # Should be in FORMAT field


def test_output_formatter_with_fragment_counts(config_vcf, temp_dir):
    """Test output with fragment counts."""
    config_vcf.output_fragment_count = True
    writer = VCFWriter(config_vcf, ["sample1"])

    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )
    variant.initialize_counts(["sample1"])
    variant.base_count["sample1"][CountType.DPF] = 5

    writer.write_variants([variant])

    # Check output includes fragment counts in FORMAT field
    with open(config_vcf.output_file) as f:
        content = f.read()
        assert "DPF" in content  # Should be in FORMAT field
