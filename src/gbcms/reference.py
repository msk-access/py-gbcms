"""Reference sequence handling."""

import logging

import pysam

logger = logging.getLogger(__name__)


class ReferenceSequence:
    """Handles reference sequence loading and access."""

    def __init__(self, fasta_file: str):
        """
        Initialize reference sequence handler.

        Args:
            fasta_file: Path to reference FASTA file (must be indexed)
        """
        self.fasta_file = fasta_file
        self.fasta: pysam.FastaFile | None = None
        self._load_reference()

    def _load_reference(self) -> None:
        """Load reference sequence using pysam."""
        logger.info(f"Loading reference sequence: {self.fasta_file}")
        try:
            self.fasta = pysam.FastaFile(self.fasta_file)
        except Exception as e:
            logger.error(f"Failed to open reference FASTA file: {e}")
            raise

    def get_base(self, chrom: str, pos: int) -> str:
        """
        Get base at specific position (0-indexed).

        Args:
            chrom: Chromosome name
            pos: 0-indexed position

        Returns:
            Base at position (uppercase)
        """
        if self.fasta is None:
            raise RuntimeError("Reference FASTA not loaded")

        try:
            return self.fasta.fetch(chrom, pos, pos + 1).upper()
        except Exception as e:
            logger.error(f"Failed to fetch base at {chrom}:{pos}: {e}")
            raise

    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        Get sequence in range (0-indexed, end exclusive).

        Args:
            chrom: Chromosome name
            start: Start position (0-indexed, inclusive)
            end: End position (0-indexed, exclusive)

        Returns:
            Sequence in range (uppercase)
        """
        if self.fasta is None:
            raise RuntimeError("Reference FASTA not loaded")

        try:
            return self.fasta.fetch(chrom, start, end).upper()
        except Exception as e:
            logger.error(f"Failed to fetch sequence at {chrom}:{start}-{end}: {e}")
            raise

    def close(self) -> None:
        """Close the FASTA file."""
        if self.fasta:
            self.fasta.close()
            self.fasta = None

    def __enter__(self) -> "ReferenceSequence":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.close()
