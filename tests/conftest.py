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
    """Create a sample FASTA file for testing."""
    fasta_file = temp_dir / "reference.fa"

    # Create a simple reference sequence
    with open(fasta_file, "w") as f:
        f.write(">chr1\n")
        f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        f.write(">chr2\n")
        f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")

    # Index the FASTA file
    pysam.faidx(str(fasta_file))

    return fasta_file


@pytest.fixture
def sample_bam(temp_dir: Path, sample_fasta: Path) -> Path:
    """Create a sample BAM file for testing."""
    bam_file = temp_dir / "sample.bam"

    # Create header
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 40},
            {"SN": "chr2", "LN": 40},
        ],
    }

    # Create BAM file
    with pysam.AlignmentFile(str(bam_file), "wb", header=header) as outf:
        # Add some test alignments
        for i in range(5):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i}"
            a.query_sequence = "ATCGATCGATCGATCGATCG"
            a.flag = 99 if i % 2 == 0 else 147  # Paired, first/second in pair
            a.reference_id = 0  # chr1
            a.reference_start = i * 2
            a.mapping_quality = 60
            a.cigarstring = "20M"  # Set CIGAR string instead of cigartuples
            a.next_reference_id = 0
            a.next_reference_start = i * 2 + 100
            a.template_length = 120
            a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIIIIIIIIIII")
            outf.write(a)

    # Index the BAM file
    pysam.index(str(bam_file))

    return bam_file


@pytest.fixture
def sample_vcf(temp_dir: Path) -> Path:
    """Create a sample VCF file for testing."""
    vcf_file = temp_dir / "variants.vcf"

    with open(vcf_file, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=40>\n")
        f.write("##contig=<ID=chr2,length=40>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("chr1\t5\t.\tA\tT\t.\tPASS\t.\n")
        f.write("chr1\t10\t.\tC\tG\t.\tPASS\t.\n")
        f.write("chr1\t15\t.\tGA\tG\t.\tPASS\t.\n")  # Deletion
        f.write("chr2\t5\t.\tG\tGC\t.\tPASS\t.\n")  # Insertion

    return vcf_file


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
        f.write("GENE1\tchr1\t5\t5\tA\tT\t\tTumor1\tNormal1\t" "10\t5\t15\t0\tMissense_Mutation\n")
        f.write("GENE2\tchr1\t10\t10\tC\tG\t\tTumor1\tNormal1\t" "8\t7\t12\t1\tMissense_Mutation\n")

    return maf_file
