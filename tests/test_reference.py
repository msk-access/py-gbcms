"""Tests for reference module."""

from pathlib import Path

import pytest

from gbcms.reference import ReferenceSequence


def test_reference_load(sample_fasta: Path):
    """Test loading reference sequence."""
    ref = ReferenceSequence(str(sample_fasta))
    assert ref.fasta is not None
    ref.close()


def test_reference_get_base(sample_fasta: Path):
    """Test getting a single base."""
    ref = ReferenceSequence(str(sample_fasta))

    # Get first base of chr1
    base = ref.get_base("chr1", 0)
    assert base == "A"

    # Get another base
    base = ref.get_base("chr1", 1)
    assert base in ["A", "T", "C", "G"]

    ref.close()


def test_reference_get_sequence(sample_fasta: Path):
    """Test getting a sequence range."""
    ref = ReferenceSequence(str(sample_fasta))

    # Get first 4 bases
    seq = ref.get_sequence("chr1", 0, 4)
    assert len(seq) == 4
    assert all(b in "ATCG" for b in seq)

    ref.close()


def test_reference_context_manager(sample_fasta: Path):
    """Test using reference as context manager."""
    with ReferenceSequence(str(sample_fasta)) as ref:
        base = ref.get_base("chr1", 0)
        assert base in ["A", "T", "C", "G"]

    # Should be closed after context
    assert ref.fasta is None


def test_reference_missing_file():
    """Test loading non-existent reference file."""
    with pytest.raises((FileNotFoundError, OSError)):
        ReferenceSequence("/nonexistent/file.fa")


def test_reference_invalid_chrom(sample_fasta: Path):
    """Test accessing invalid chromosome."""
    ref = ReferenceSequence(str(sample_fasta))

    with pytest.raises(KeyError):
        ref.get_base("chrNonexistent", 0)

    ref.close()


def test_reference_invalid_position(sample_fasta: Path):
    """Test accessing invalid position."""
    ref = ReferenceSequence(str(sample_fasta))

    # Position beyond chromosome length should fail or return empty
    # Note: pysam behavior may vary, so we test that it doesn't crash
    try:
        result = ref.get_base("chr1", 10000)
        # If it doesn't raise an exception, result should be empty or None
        assert result == "" or result is None
    except (IndexError, ValueError):
        # Some versions of pysam raise IndexError or ValueError
        pass

    ref.close()
