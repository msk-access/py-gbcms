"""Tests for variant module."""

from pathlib import Path

import pytest

from gbcms.variant import VariantEntry, VariantLoader


def test_variant_entry_creation():
    """Test creating a VariantEntry."""
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )

    assert variant.chrom == "chr1"
    assert variant.pos == 100
    assert variant.ref == "A"
    assert variant.alt == "T"
    assert variant.snp is True
    assert variant.insertion is False
    assert variant.deletion is False


def test_variant_entry_key():
    """Test variant key generation."""
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
    )

    key = variant.get_variant_key()
    assert key == ("chr1", 100, "A", "T")


def test_variant_entry_sorting():
    """Test variant sorting."""
    v1 = VariantEntry(chrom="chr1", pos=100, end_pos=100, ref="A", alt="T")
    v2 = VariantEntry(chrom="chr1", pos=200, end_pos=200, ref="C", alt="G")
    v3 = VariantEntry(chrom="chr2", pos=50, end_pos=50, ref="G", alt="A")

    variants = [v2, v3, v1]
    variants.sort()

    assert variants[0] == v1
    assert variants[1] == v2
    assert variants[2] == v3


def test_variant_loader_vcf(sample_vcf: Path):
    """Test loading variants from VCF."""
    loader = VariantLoader()
    variants = loader.load_vcf(str(sample_vcf))

    assert len(variants) == 4

    # Check first variant (SNP)
    assert variants[0].chrom == "chr1"
    assert variants[0].pos == 4  # 0-indexed
    assert variants[0].ref == "A"
    assert variants[0].alt == "T"
    assert variants[0].snp is True

    # Check deletion
    assert variants[2].deletion is True
    assert variants[2].ref == "GA"
    assert variants[2].alt == "G"

    # Check insertion
    assert variants[3].insertion is True
    assert variants[3].ref == "G"
    assert variants[3].alt == "GC"


def test_variant_loader_maf(sample_maf: Path):
    """Test loading variants from MAF."""

    # Mock reference getter
    def mock_ref_getter(chrom: str, pos: int) -> str:
        return "A"

    loader = VariantLoader(reference_getter=mock_ref_getter)
    variants = loader.load_maf(str(sample_maf))

    assert len(variants) == 2

    # Check first variant
    assert variants[0].chrom == "chr1"
    assert variants[0].pos == 4  # 0-indexed
    assert variants[0].ref == "A"
    assert variants[0].alt == "T"
    assert variants[0].gene == "GENE1"
    assert variants[0].tumor_sample == "Tumor1"
    assert variants[0].normal_sample == "Normal1"
    assert variants[0].t_ref_count == 10
    assert variants[0].t_alt_count == 5


def test_variant_loader_maf_missing_columns(temp_dir: Path):
    """Test loading MAF with missing required columns."""
    maf_file = temp_dir / "bad.maf"
    with open(maf_file, "w") as f:
        f.write("Hugo_Symbol\tChromosome\n")
        f.write("GENE1\tchr1\n")

    loader = VariantLoader()
    with pytest.raises(ValueError, match="Incorrect MAF file header"):
        loader.load_maf(str(maf_file))


def test_variant_entry_initialize_counts():
    """Test initializing counts for samples."""
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
    )

    samples = ["sample1", "sample2"]
    variant.initialize_counts(samples)

    assert "sample1" in variant.base_count
    assert "sample2" in variant.base_count
    assert len(variant.base_count["sample1"]) == 9  # 9 count types


def test_variant_entry_get_count():
    """Test getting counts from variant."""
    variant = VariantEntry(
        chrom="chr1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
    )

    variant.initialize_counts(["sample1"])

    # Default should be 0
    from gbcms.config import CountType

    assert variant.get_count("sample1", CountType.DP) == 0.0

    # Set a value
    variant.base_count["sample1"][CountType.DP] = 10.0
    assert variant.get_count("sample1", CountType.DP) == 10.0

    # Non-existent sample should return 0
    assert variant.get_count("nonexistent", CountType.DP) == 0.0
