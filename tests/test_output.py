"""Tests for output module."""

from pathlib import Path

import pytest

from gbcms.config import Config, CountType
from gbcms.output import OutputFormatter
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
        output_maf=True,
    )


def test_output_formatter_vcf(config_vcf, temp_dir):
    """Test VCF output formatting."""
    formatter = OutputFormatter(config_vcf, ["sample1", "sample2"])

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

    formatter.write_vcf_output([variant])

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


def test_output_formatter_maf(config_maf, temp_dir):
    """Test MAF output formatting."""
    formatter = OutputFormatter(config_maf, ["Tumor1", "Normal1"])

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
        maf_pos=100,
        maf_end_pos=100,
        maf_ref="A",
        maf_alt="T",
        maf_line="GENE1\tchr1\t101\t101\tA\tT\t\tTumor1\tNormal1\tMissense_Mutation\t10\t6\t4\t15\t15\t0\t0\t0\t0\t0\t0\t0",  # Mock MAF line
    )
    variant.initialize_counts(["Tumor1", "Normal1"])
    variant.base_count["Tumor1"][CountType.DP] = 10
    variant.base_count["Tumor1"][CountType.RD] = 6
    variant.base_count["Tumor1"][CountType.AD] = 4
    variant.base_count["Normal1"][CountType.DP] = 15
    variant.base_count["Normal1"][CountType.RD] = 15
    variant.base_count["Normal1"][CountType.AD] = 0

    formatter.write_maf_output([variant])

    # Check output file exists
    assert Path(config_maf.output_file).exists()

    # Read and verify content
    with open(config_maf.output_file) as f:
        lines = f.readlines()
        assert len(lines) == 2  # Header + 1 variant
        assert "Hugo_Symbol" in lines[0]
        assert "GENE1" in lines[1]


def test_output_formatter_fillout(config_maf, temp_dir):
    """Test fillout output formatting."""
    config_maf.output_maf = False
    formatter = OutputFormatter(config_maf, ["Tumor1", "Normal1"])

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
        maf_pos=100,
        maf_end_pos=100,
        maf_ref="A",
        maf_alt="T",
    )
    variant.initialize_counts(["Tumor1", "Normal1"])

    formatter.write_fillout_output([variant])

    # Check output file exists
    assert Path(config_maf.output_file).exists()


def test_output_formatter_with_strand_counts(config_vcf, temp_dir):
    """Test output with strand counts."""
    config_vcf.output_positive_count = True
    formatter = OutputFormatter(config_vcf, ["sample1"])

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
    variant.base_count["sample1"][CountType.DPP] = 6

    formatter.write_vcf_output([variant])

    # Check output includes strand counts in FORMAT field
    with open(config_vcf.output_file) as f:
        content = f.read()
        assert "DPP" in content  # Should be in FORMAT field


def test_output_formatter_with_fragment_counts(config_vcf, temp_dir):
    """Test output with fragment counts."""
    config_vcf.output_fragment_count = True
    formatter = OutputFormatter(config_vcf, ["sample1"])

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

    formatter.write_vcf_output([variant])

    # Check output includes fragment counts in FORMAT field
    with open(config_vcf.output_file) as f:
        content = f.read()
        assert "DPF" in content  # Should be in FORMAT field
