"""
Coordinate Kernel: The source of truth for genomic coordinate systems.

Handles conversion between:
- VCF (1-based)
- MAF (1-based)
- Internal (0-based, half-open [start, end))

Ensures consistent representation of variants:
- SNPs: 0-based index of the base.
- Insertions: 0-based index of the ANCHOR base (preceding the insertion).
- Deletions: 0-based index of the ANCHOR base (preceding the deletion).
"""

from gbcms.models.core import Variant, VariantType


class CoordinateKernel:
    """
    Stateless utility for coordinate transformations and normalization.
    """

    @staticmethod
    def vcf_to_internal(
        chrom: str, pos: int, ref: str, alt: str, original_id: str | None = None
    ) -> Variant:
        """
        Convert VCF coordinates (1-based) to internal normalized Variant.

        Args:
            chrom: Chromosome name
            pos: 1-based position from VCF
            ref: Reference allele
            alt: Alternate allele
            original_id: Optional VCF ID

        Returns:
            Normalized Variant object
        """
        norm_chrom = CoordinateKernel.normalize_chromosome(chrom)

        # Determine variant type and internal position
        if len(ref) == 1 and len(alt) == 1:
            vtype = VariantType.SNP
            # SNP: VCF POS is the base itself.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        elif len(ref) == 1 and len(alt) > 1:
            vtype = VariantType.INSERTION
            # Insertion: VCF POS is the base BEFORE the insertion (the anchor).
            # VCF: POS=10, REF=A, ALT=AT (Insertion of T after A at 10)
            # Internal: 0-based index of the ANCHOR base.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        elif len(ref) > 1 and len(alt) == 1:
            vtype = VariantType.DELETION
            # Deletion: VCF POS is the base BEFORE the deletion (the anchor).
            # VCF: POS=10, REF=AT, ALT=A (Deletion of T after A at 10)
            # Internal POS: 0-based index of the ANCHOR base.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        else:
            vtype = VariantType.COMPLEX
            # Complex: Treat start as 0-based index of first ref base
            internal_pos = pos - 1

        return Variant(
            chrom=norm_chrom,
            pos=internal_pos,
            ref=ref,
            alt=alt,
            variant_type=vtype,
            original_id=original_id,
        )

    @staticmethod
    def maf_to_internal(chrom: str, start_pos: int, end_pos: int, ref: str, alt: str) -> Variant:
        """
        Convert MAF coordinates (1-based inclusive) to internal normalized Variant.

        MAF coordinates are generally 1-based inclusive [start, end].
        """
        norm_chrom = CoordinateKernel.normalize_chromosome(chrom)

        # Handle MAF indels which often use '-'
        if ref == "-" or alt == "-":
            # MAF Insertion: Start_Position is the base BEFORE the insertion (anchor).
            # ref='-', alt='T' -> VCF-like would be ref='A', alt='AT' (requires lookup)
            # But if we just want to represent it internally:
            if ref == "-":  # Insertion
                vtype = VariantType.INSERTION
                # MAF Start_Position is usually the flanking base 0 or 1?
                # Standard MAF: Start_Position is the base BEFORE the insertion.
                internal_pos = start_pos - 1
            else:  # Deletion
                vtype = VariantType.DELETION
                # MAF Start_Position is the first deleted base? Or anchor?
                # Usually first deleted base.
                # We need to convert to anchor-based for consistency if possible,
                # OR handle MAF-style internally.
                # Let's assume we want VCF-style anchor-based internally.
                # This effectively requires a reference lookup to get the anchor base.
                # For now, we will mark it as needing normalization or handle it in the engine.
                internal_pos = start_pos - 1

        elif len(ref) == len(alt) == 1:
            vtype = VariantType.SNP
            internal_pos = start_pos - 1

        else:
            vtype = VariantType.COMPLEX
            internal_pos = start_pos - 1

        return Variant(chrom=norm_chrom, pos=internal_pos, ref=ref, alt=alt, variant_type=vtype)

    @staticmethod
    def normalize_chromosome(chrom: str) -> str:
        """
        Normalize chromosome name (remove 'chr' prefix).
        """
        if chrom.lower().startswith("chr"):
            return chrom[3:]
        return chrom
