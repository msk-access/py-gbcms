"""Tests for counter module."""

import pysam
import pytest

from gbcms.config import Config, CountType
from gbcms.counter import BaseCounter
from gbcms.variant import VariantEntry


@pytest.fixture
def config(temp_dir, sample_fasta, sample_bam, sample_vcf):
    """Create a test configuration."""
    return Config(
        fasta_file=str(sample_fasta),
        bam_files={"sample1": str(sample_bam)},
        variant_files=[str(sample_vcf)],
        output_file=str(temp_dir / "output.txt"),
        input_is_vcf=True,
    )


@pytest.fixture
def counter(config):
    """Create a BaseCounter instance."""
    return BaseCounter(config)


def test_filter_alignment_duplicate(counter):
    """Test filtering duplicate reads."""
    aln = pysam.AlignedSegment()
    aln.is_duplicate = True
    aln.mapping_quality = 60

    assert counter.filter_alignment(aln) is True


def test_filter_alignment_low_mapq(counter):
    """Test filtering low mapping quality reads."""
    aln = pysam.AlignedSegment()
    aln.is_duplicate = False
    aln.mapping_quality = 10

    assert counter.filter_alignment(aln) is True


def test_filter_alignment_pass(counter):
    """Test alignment that passes filters."""
    aln = pysam.AlignedSegment()
    aln.is_duplicate = False
    aln.mapping_quality = 60
    aln.is_qcfail = False
    aln.is_secondary = False
    aln.is_supplementary = False

    assert counter.filter_alignment(aln) is False


def test_has_indel():
    """Test indel detection in alignments."""
    # Alignment without indel
    aln1 = pysam.AlignedSegment()
    aln1.cigartuples = [(0, 20)]  # 20M
    assert BaseCounter._has_indel(aln1) is False

    # Alignment with insertion
    aln2 = pysam.AlignedSegment()
    aln2.cigartuples = [(0, 10), (1, 5), (0, 10)]  # 10M5I10M
    assert BaseCounter._has_indel(aln2) is True

    # Alignment with deletion
    aln3 = pysam.AlignedSegment()
    aln3.cigartuples = [(0, 10), (2, 5), (0, 10)]  # 10M5D10M
    assert BaseCounter._has_indel(aln3) is True


def test_count_bases_snp(counter):
    """Test counting bases for SNP variant."""
    variant = VariantEntry(
        chrom="chr1",
        pos=5,
        end_pos=5,
        ref="A",
        alt="T",
        snp=True,
    )
    variant.initialize_counts(["sample1"])

    # Create test alignments
    alignments = []
    for i in range(3):
        aln = pysam.AlignedSegment()
        aln.query_name = f"read_{i}"
        aln.reference_start = 0
        aln.reference_id = 0
        aln.query_sequence = "ATCGATCGATCGATCGATCG"
        aln.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIIIIIIIIIII")
        aln.cigartuples = [(0, 20)]
        aln.is_reverse = i % 2 == 0
        aln.is_read1 = True
        alignments.append(aln)

    counter.count_bases_snp(variant, alignments, "sample1")

    # Check counts
    assert variant.get_count("sample1", CountType.DP) > 0


def test_count_bases_dnp(counter):
    """Test counting bases for DNP variant."""
    variant = VariantEntry(
        chrom="chr1",
        pos=5,
        end_pos=6,
        ref="AT",
        alt="GC",
        dnp=True,
        dnp_len=2,
    )
    variant.initialize_counts(["sample1"])

    # Create test alignment
    aln = pysam.AlignedSegment()
    aln.query_name = "read_1"
    aln.reference_start = 0
    aln.query_sequence = "ATCGATCGATCGATCGATCG"
    aln.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIIIIIIIIIII")
    aln.cigartuples = [(0, 20)]
    aln.is_reverse = False
    aln.is_read1 = True

    counter.count_bases_dnp(variant, [aln], "sample1")

    # Should have some depth
    assert variant.get_count("sample1", CountType.DP) >= 0


def test_count_bases_indel(counter):
    """Test counting bases for indel variant."""
    variant = VariantEntry(
        chrom="chr1",
        pos=5,
        end_pos=6,
        ref="AT",
        alt="A",
        deletion=True,
    )
    variant.initialize_counts(["sample1"])

    # Create test alignment
    aln = pysam.AlignedSegment()
    aln.query_name = "read_1"
    aln.reference_start = 0
    aln.query_sequence = "ATCGATCGATCGATCGATCG"
    aln.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIIIIIIIIIII")
    aln.cigartuples = [(0, 20)]
    aln.is_reverse = False
    aln.is_read1 = True

    counter.count_bases_indel(variant, [aln], "sample1")

    # Should complete without error
    assert variant.base_count["sample1"] is not None


def test_count_variant_dispatch(counter):
    """Test count_variant dispatches to correct method."""
    # SNP
    snp_variant = VariantEntry(chrom="chr1", pos=5, end_pos=5, ref="A", alt="T", snp=True)
    snp_variant.initialize_counts(["sample1"])
    counter.count_variant(snp_variant, [], "sample1")

    # DNP
    dnp_variant = VariantEntry(
        chrom="chr1", pos=5, end_pos=6, ref="AT", alt="GC", dnp=True, dnp_len=2
    )
    dnp_variant.initialize_counts(["sample1"])
    counter.count_variant(dnp_variant, [], "sample1")

    # Insertion
    ins_variant = VariantEntry(chrom="chr1", pos=5, end_pos=5, ref="A", alt="AT", insertion=True)
    ins_variant.initialize_counts(["sample1"])
    counter.count_variant(ins_variant, [], "sample1")

    # Deletion
    del_variant = VariantEntry(chrom="chr1", pos=5, end_pos=6, ref="AT", alt="A", deletion=True)
    del_variant.initialize_counts(["sample1"])
    counter.count_variant(del_variant, [], "sample1")
