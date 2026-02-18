"""
Input Adapters: Handling VCF and MAF inputs.

This module provides classes to read variants from VCF and MAF files,
converting them into the internal normalized representation using CoordinateKernel.
"""

import csv
import logging
import warnings
from collections.abc import Iterator
from pathlib import Path

import pysam
from pydantic import ValidationError

from ..core.kernel import CoordinateKernel
from ..models.core import Variant

logger = logging.getLogger(__name__)

__all__ = ["VariantReader", "VcfReader", "MafReader", "ReferenceChecker"]


class VariantReader:
    """Abstract base class for variant readers."""

    def __iter__(self) -> Iterator[Variant]:
        raise NotImplementedError


class VcfReader(VariantReader):
    """Reads variants from a VCF file."""

    def __init__(self, path: Path):
        self.path = path
        self._vcf = pysam.VariantFile(str(path))

    def __iter__(self) -> Iterator[Variant]:
        for record in self._vcf:
            # VCF coordinates are 1-based
            # pysam converts them to 0-based automatically?
            # pysam.VariantFile returns 0-based pos (start)
            # BUT CoordinateKernel.vcf_to_internal expects 1-based VCF POS.
            # Let's check pysam documentation or behavior.
            # pysam record.pos is 0-based. record.start is 0-based.
            # The VCF file itself has 1-based POS.
            # If we use record.pos + 1, we get the VCF POS.

            # Handle multiple ALTs
            for alt in record.alts or []:
                # VCF POS is record.pos (1-based) or record.start + 1
                if not record.ref:
                    continue  # Skip if no REF

                yield CoordinateKernel.vcf_to_internal(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alt=alt,
                    original_id=record.id,
                )

    def close(self):
        self._vcf.close()


class MafReader(VariantReader):
    """Reads variants from a MAF file.

    Normalizes MAF-style indels and complex variants to VCF-style
    representation by fetching anchor bases from the reference genome.

    Args:
        path: Path to the MAF file.
        fasta_path: Optional path to reference FASTA for indel normalization.
    """

    def __init__(self, path: Path, fasta_path: Path | None = None):
        self.path = path
        self.fasta = pysam.FastaFile(str(fasta_path)) if fasta_path else None

    def _fetch_anchor_base(self, chrom: str, pos_0based: int) -> str | None:
        """Fetch a single reference base, trying normalized then original chrom name.

        Args:
            chrom: Chromosome name (unnormalized).
            pos_0based: 0-based genomic position.

        Returns:
            Uppercase base character, or None if fetch fails.
        """
        if not self.fasta:
            return None
        norm_chrom = CoordinateKernel.normalize_chromosome(chrom)
        for name in (norm_chrom, chrom):
            try:
                return self.fasta.fetch(name, pos_0based, pos_0based + 1).upper()
            except (KeyError, ValueError):
                continue
        logger.warning(
            "FASTA fetch failed at %s:%d (tried '%s' and '%s')",
            chrom,
            pos_0based,
            norm_chrom,
            chrom,
        )
        return None

    def __iter__(self) -> Iterator[Variant]:
        with open(self.path) as f:
            # Skip comments
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith("#"):
                    f.seek(pos)
                    break

            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                try:
                    chrom = row["Chromosome"]
                    start_pos = int(row["Start_Position"])
                    ref = row["Reference_Allele"]
                    alt = row["Tumor_Seq_Allele2"]  # Standard MAF alt column

                    # Normalize indels and complex variants if FASTA is available.
                    # Condition: normalize anything that isn't a simple SNP.
                    is_snp = len(ref) == 1 and len(alt) == 1 and ref != "-" and alt != "-"
                    if self.fasta and not is_snp:
                        if ref == "-":  # Insertion
                            # MAF Start_Position is the base BEFORE the insertion
                            anchor_pos_1based = start_pos
                            anchor_base = self._fetch_anchor_base(chrom, anchor_pos_1based - 1)
                            if anchor_base is None:
                                continue

                            # VCF: POS=anchor, REF=anchor_base, ALT=anchor_base+ins
                            vcf_pos = anchor_pos_1based
                            vcf_ref = anchor_base
                            vcf_alt = anchor_base + alt

                        elif alt == "-":  # Deletion
                            # MAF Start_Position is the FIRST DELETED base
                            anchor_pos_1based = start_pos - 1
                            anchor_base = self._fetch_anchor_base(chrom, anchor_pos_1based - 1)
                            if anchor_base is None:
                                continue

                            # VCF: POS=anchor, REF=anchor_base+del, ALT=anchor_base
                            vcf_pos = anchor_pos_1based
                            vcf_ref = anchor_base + ref
                            vcf_alt = anchor_base

                        else:  # Complex: both non-dash, different lengths
                            # MAF Start_Position is first changed base;
                            # anchor is the base before that.
                            anchor_pos_1based = start_pos - 1
                            anchor_base = self._fetch_anchor_base(chrom, anchor_pos_1based - 1)
                            if anchor_base is None:
                                continue

                            vcf_pos = anchor_pos_1based
                            vcf_ref = anchor_base + ref
                            vcf_alt = anchor_base + alt
                            logger.debug(
                                "Complex variant normalized: %s:%d %s>%s " "→ VCF %s:%d %s>%s",
                                chrom,
                                start_pos,
                                ref,
                                alt,
                                chrom,
                                vcf_pos,
                                vcf_ref,
                                vcf_alt,
                            )

                        yield CoordinateKernel.vcf_to_internal(
                            chrom=chrom, pos=vcf_pos, ref=vcf_ref, alt=vcf_alt
                        ).model_copy(update={"metadata": row})
                    else:
                        # Fallback to old behavior or direct mapping for SNPs
                        # For SNPs, MAF Start_Position == VCF POS
                        if len(ref) == len(alt) == 1 and ref != "-" and alt != "-":
                            yield CoordinateKernel.vcf_to_internal(
                                chrom=chrom, pos=start_pos, ref=ref, alt=alt
                            ).model_copy(update={"metadata": row})
                        else:
                            # Fallback for complex/unhandled without FASTA.
                            # This path may produce incorrect indel coordinates
                            # because maf_to_internal cannot properly resolve
                            # anchor bases without a reference genome.
                            if ref == "-" or alt == "-":
                                logger.warning(
                                    "MAF indel at %s:%s uses maf_to_internal() fallback "
                                    "(no --fasta provided). Indel coordinates may be incorrect. "
                                    "Provide --fasta for accurate indel normalization.",
                                    chrom,
                                    start_pos,
                                )
                                warnings.warn(
                                    "maf_to_internal() is deprecated for indels. "
                                    "Use --fasta for proper VCF-style normalization.",
                                    DeprecationWarning,
                                    stacklevel=2,
                                )
                            yield CoordinateKernel.maf_to_internal(
                                chrom=chrom,
                                start_pos=start_pos,
                                end_pos=int(row["End_Position"]),
                                ref=ref,
                                alt=alt,
                            ).model_copy(update={"metadata": row})

                except (KeyError, ValueError, ValidationError):
                    # Log warning or skip malformed lines
                    continue

    def close(self):
        if self.fasta:
            self.fasta.close()


class ReferenceChecker:
    """
    Utility to check variants against a reference FASTA.
    Ensures that the REF allele matches the genome.
    """

    def __init__(self, fasta_path: Path):
        self.fasta = pysam.FastaFile(str(fasta_path))

    def validate(self, variant: Variant) -> bool:
        """
        Check if variant REF matches reference genome.
        """
        # Variant pos is 0-based.
        # Fetch sequence of length REF
        try:
            # Try normalized and potentially 'chr' prefixed chromosome names
            chrom = variant.chrom
            # chrom is already normalized (e.g. "1") by CoordinateKernel

            ref_seq = None
            try:
                ref_seq = self.fasta.fetch(chrom, variant.pos, variant.pos + len(variant.ref))
            except (ValueError, KeyError):
                try:
                    # Try adding 'chr' prefix
                    ref_seq = self.fasta.fetch(
                        f"chr{chrom}", variant.pos, variant.pos + len(variant.ref)
                    )
                except (ValueError, KeyError) as e:
                    logger.debug("Failed to fetch %s and chr%s: %s", chrom, chrom, e)
                    return False

            if ref_seq is None:
                return False

            return ref_seq.upper() == variant.ref.upper()

        except Exception:
            return False

    def close(self):
        self.fasta.close()
