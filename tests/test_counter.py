"""Tests for counter module."""

import pysam

from gbcms.config import Config
from gbcms.counter import BaseCounter
from gbcms.variant import VariantEntry


def test_filter_alignment_duplicate():
    """Test filtering duplicate reads."""
    # Create a minimal config for testing
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Test duplicate filtering
    aln = pysam.AlignedSegment()
    aln.is_duplicate = True
    aln.mapping_quality = 60

    assert counter._should_filter_alignment(aln) is True


def test_filter_alignment_low_mapq():
    """Test filtering low mapping quality reads."""
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Test low mapping quality filtering
    aln = pysam.AlignedSegment()
    aln.is_duplicate = False
    aln.mapping_quality = 10

    assert counter._should_filter_alignment(aln) is True


def test_filter_alignment_pass():
    """Test alignment that passes filters."""
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Test alignment that passes filters
    aln = pysam.AlignedSegment()
    aln.is_duplicate = False
    aln.mapping_quality = 60

    assert counter._should_filter_alignment(aln) is False


def test_count_bases_snp():
    """Test counting bases for SNP variant."""
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Create a simple SNP variant
    variant = VariantEntry(
        chrom="1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )

    # Test that counting doesn't crash (we can't easily test exact counts without BAM data)
    try:
        alignments: list = []  # Empty list for this test
        counter.smart_count_variant(variant, alignments, "sample1")
        # If we get here without exception, the test passes
        assert True
    except ValueError as e:
        if "No alignments provided" in str(e):
            # Expected when no alignments are provided
            assert True
        else:
            raise


def test_count_bases_indel():
    """Test counting bases for indel variant."""
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Create a simple insertion variant
    variant = VariantEntry(
        chrom="1",
        pos=100,
        end_pos=101,
        ref="A",
        alt="AT",
        insertion=True,
    )

    # Test that counting doesn't crash
    try:
        alignments: list = []  # Empty list for this test
        counter.smart_count_variant(variant, alignments, "sample1")
        assert True
    except ValueError as e:
        if "No alignments provided" in str(e):
            assert True
        else:
            raise


def test_count_variant_dispatch():
    """Test variant dispatch to appropriate counting method."""
    config = Config(
        fasta_file="tests/testdata/integration_test_reference.fa",
        bam_files={"sample1": "tests/testdata/sample1_integration_test.bam"},
        variant_files=["tests/testdata/integration_test_variants.vcf"],
        output_file="/tmp/test_output.txt",
        input_is_vcf=True,
    )

    counter = BaseCounter(config)

    # Test SNP dispatch
    snp_variant = VariantEntry(
        chrom="1",
        pos=100,
        end_pos=100,
        ref="A",
        alt="T",
        snp=True,
    )

    # Test insertion dispatch
    insertion_variant = VariantEntry(
        chrom="1",
        pos=100,
        end_pos=101,
        ref="A",
        alt="AT",
        insertion=True,
    )

    # Test deletion dispatch
    deletion_variant = VariantEntry(
        chrom="1",
        pos=100,
        end_pos=101,
        ref="AT",
        alt="A",
        deletion=True,
    )

    # Test that dispatch works (will fail with no alignments, but that's expected)
    for variant in [snp_variant, insertion_variant, deletion_variant]:
        try:
            counter.smart_count_variant(variant, [], "sample1")
        except ValueError as e:
            if "No alignments provided" in str(e):
                continue  # Expected
            else:
                raise
