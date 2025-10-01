"""Tests for reference module."""

import pytest
from pathlib import Path

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
    with pytest.raises(Exception):
        ReferenceSequence("/nonexistent/file.fa")


def test_reference_invalid_chrom(sample_fasta: Path):
    """Test accessing invalid chromosome."""
    ref = ReferenceSequence(str(sample_fasta))
    
    with pytest.raises(Exception):
        ref.get_base("chrNonexistent", 0)
    
    ref.close()


def test_reference_invalid_position(sample_fasta: Path):
    """Test accessing invalid position."""
    ref = ReferenceSequence(str(sample_fasta))
    
    # Position beyond chromosome length should fail or return empty
    with pytest.raises(Exception):
        ref.get_base("chr1", 10000)
    
    ref.close()
