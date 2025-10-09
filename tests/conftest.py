"""Pytest configuration and fixtures."""

import sys
import tempfile
from collections.abc import Generator
from pathlib import Path

import pysam
import pytest

# Add src directory to path so tests use local code, not installed package
project_root = Path(__file__).parent.parent
src_dir = project_root / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fasta(temp_dir: Path) -> Path:
    """Use real genomic data for testing."""
    test_data_dir = Path(__file__).parent / "testdata"
    return test_data_dir / "integration_test_reference.fa"


@pytest.fixture
def sample_bam(temp_dir: Path, sample_fasta: Path) -> Path:
    """Use real BAM data for testing."""
    test_data_dir = Path(__file__).parent / "testdata"
    return test_data_dir / "sample1_integration_test.bam"


@pytest.fixture
def sample_vcf(temp_dir: Path) -> Path:
    """Use real VCF data for testing."""
    test_data_dir = Path(__file__).parent / "testdata"
    return test_data_dir / "integration_test_variants.vcf"


@pytest.fixture
def sample_maf(temp_dir: Path) -> Path:
    """Create a sample MAF file for testing."""
    maf_file = temp_dir / "variants.maf"

    with open(maf_file, "w") as f:
        # Header
        f.write(
            "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
            "Reference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\t"
            "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\t"
            "t_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\t"
            "Variant_Classification\n"
        )
        # Variants
        f.write("GENE1\tchr1\t5\t5\tA\tA\tT\tTumor1\tNormal1\t" "10\t5\t15\t0\tMissense_Mutation\n")
        f.write("GENE2\tchr1\t10\t10\tC\tC\tG\tTumor1\tNormal1\t" "8\t7\t12\t1\tMissense_Mutation\n")

    return maf_file


@pytest.fixture
def multi_sample_bams(temp_dir: Path) -> dict:
    """Provide multiple BAM files for multi-sample testing."""
    test_data_dir = Path(__file__).parent / "testdata"
    return {
        "sample1": test_data_dir / "sample1_integration_test.bam",
        "sample2": test_data_dir / "sample2_integration_test.bam",
    }
