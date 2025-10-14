"""Variant loading and representation."""

import logging
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .config import CountType

logger = logging.getLogger(__name__)

# Try to import cyvcf2 for fast VCF parsing
try:
    from cyvcf2 import VCF

    HAS_CYVCF2 = True
    logger.debug("cyvcf2 available - using fast VCF parsing")
except ImportError:
    HAS_CYVCF2 = False
    logger.debug("cyvcf2 not available - using pure Python VCF parsing")


def _is_standard_chromosome(chrom: str) -> bool:
    """
    Check if chromosome should be normalized.

    Only standard chromosomes (1-22, X, Y, M, MT) are normalized.
    Alternative contigs (KI270728.1, etc.) are left unchanged.

    Args:
        chrom: Chromosome name

    Returns:
        True if chromosome should be normalized
    """
    if not chrom:
        return False

    # Standard chromosomes that can be normalized (check after removing chr prefix)
    standard_chromosomes = {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "M",
        "MT",
    }

    # Case-insensitive removal of chr prefix and check
    base_chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE).upper()
    return base_chrom in standard_chromosomes


def normalize_chromosome_name_to_format(chrom: str, target_format: str) -> str:
    """
    Normalize chromosome name to match target format (bidirectional).

    This function can add or strip "chr" prefix based on the target format,
    and only applies to standard chromosomes (1-22, X, Y, M, MT).

    Args:
        chrom: Original chromosome name
        target_format: "chr_prefix" or "no_prefix"

    Returns:
        Normalized chromosome name
    """
    if not chrom:
        return chrom

    # Only normalize standard chromosomes
    if not _is_standard_chromosome(chrom):
        return chrom

    if target_format == "chr_prefix":
        # Add chr prefix if not present
        return f"chr{chrom}" if not chrom.startswith("chr") else chrom
    else:  # no_prefix
        # Strip chr prefix if present
        return chrom[3:] if chrom.startswith("chr") else chrom


def normalize_chromosome_name(chrom: str) -> str:
    """
    Normalize chromosome names to a consistent format.

    Handles common variations:
    - "chr1" <-> "1"
    - "chrX" <-> "X"
    - "chrY" <-> "Y"
    - "chrM" <-> "MT"

    Args:
        chrom: Original chromosome name

    Returns:
        Normalized chromosome name
    """
    if not chrom:
        return chrom

    # Remove chr prefix if present
    if chrom.startswith("chr"):
        chrom = chrom[3:]

    # Handle mitochondrial chromosome variations
    if chrom in ("M", "MT", "m", "mt"):
        return "MT"

    # Handle sex chromosomes
    if chrom in ("X", "Y"):
        return chrom

    # Validate numeric chromosomes
    try:
        int(chrom)  # Should be 1-22 or numeric string
        return chrom
    except ValueError:
        logger.warning(f"Unexpected chromosome name format: {chrom}")
        return chrom


@dataclass
class VariantEntry:
    """Represents a variant with its counts across samples."""

    chrom: str
    bam_pos: int = 0  # 0-indexed for PySAM BAM access
    bam_end_pos: int = 0  # 0-indexed end position for PySAM
    pos: int = 0  # Alias for bam_pos for backward compatibility
    end_pos: int = 0  # Alias for bam_end_pos for backward compatibility
    ref: str = ""
    alt: str = ""
    snp: bool = False

    insertion: bool = False
    deletion: bool = False
    tumor_sample: str = ""
    normal_sample: str = ""
    gene: str = ""
    effect: str = ""
    t_ref_count: int = 0
    t_alt_count: int = 0
    n_ref_count: int = 0
    n_alt_count: int = 0
    original_pos: int = 0  # Original 1-indexed coordinates from input
    original_end_pos: int = 0  # Original end position from input
    maf_pos: int = 0  # MAF-specific position
    maf_end_pos: int = 0  # MAF-specific end position
    maf_ref: str = ""
    maf_alt: str = ""
    caller: str = ""
    base_count: dict[str, np.ndarray] = field(default_factory=dict)
    duplicate_variant_ptr: Optional["VariantEntry"] = None
    maf_line: str = ""  # Store original MAF line for output

    def __post_init__(self):
        """Synchronize coordinate fields for backward compatibility."""
        # Sync pos <-> bam_pos
        if self.pos and not self.bam_pos:
            object.__setattr__(self, 'bam_pos', self.pos)
        elif self.bam_pos and not self.pos:
            object.__setattr__(self, 'pos', self.bam_pos)
            
        # Sync end_pos <-> bam_end_pos  
        if self.end_pos and not self.bam_end_pos:
            object.__setattr__(self, 'bam_end_pos', self.end_pos)
        elif self.bam_end_pos and not self.end_pos:
            object.__setattr__(self, 'end_pos', self.bam_end_pos)



    def get_variant_key(self) -> tuple[str, int, str, str]:
        """Return unique key for variant identification."""
        return (self.chrom, self.bam_pos, self.ref, self.alt)

    def initialize_counts(self, sample_names: list[str]) -> None:
        """Initialize count arrays for all samples."""
        for sample in sample_names:
            if sample not in self.base_count:
                self.base_count[sample] = np.zeros(len(CountType), dtype=np.float32)

    def get_count(self, sample: str, count_type: CountType) -> float:
        """Get count for specific sample and type."""
        if sample not in self.base_count:
            return 0.0
        return float(self.base_count[sample][count_type])

    def __lt__(self, other: "VariantEntry") -> bool:
        """Compare variants for sorting."""
        if self.chrom != other.chrom:
            return self._chrom_sort_key() < other._chrom_sort_key()
        return self.bam_pos < other.bam_pos

    def _chrom_sort_key(self) -> tuple:
        """Generate sort key for chromosome."""
        # Chromosomes are now already normalized to target format during loading
        # No need to remove "chr" prefix as it's handled during normalization
        chrom = self.chrom
        try:
            return (0, int(chrom))
        except ValueError:
            if chrom == "X":
                return (1, 0)
            elif chrom == "Y":
                return (1, 1)
            elif chrom == "M" or chrom == "MT":
                return (1, 2)
            else:
                return (2, chrom)


class VariantLoader:
    """Loads variants from VCF or MAF files."""

    def __init__(
        self, reference_getter=None, chromosome_normalization_map=None, target_format=None
    ):
        """
        Initialize variant loader.

        Args:
            reference_getter: Callable that takes (chrom, pos) and returns base
            chromosome_normalization_map: Dict mapping original -> normalized chromosome names
            target_format: Target chromosome format ("chr_prefix" or "no_prefix")
        """
        self.reference_getter = reference_getter
        self.chromosome_normalization_map: dict[str, str] = chromosome_normalization_map or {}
        self.target_format = target_format

    def _normalize_chromosome(self, chrom: str) -> str:
        """Normalize chromosome name using the provided normalization map."""
        return self.chromosome_normalization_map.get(chrom, chrom)

    def load_vcf(self, vcf_file: str) -> list[VariantEntry]:
        """
        Load variants from VCF file.

        Uses cyvcf2 for fast parsing if available, otherwise falls back to pure Python.

        Args:
            vcf_file: Path to VCF file (can be .vcf or .vcf.gz)

        Returns:
            List of VariantEntry objects
        """
        if HAS_CYVCF2:
            return self._load_vcf_cyvcf2(vcf_file)
        else:
            return self._load_vcf_python(vcf_file)

    def _load_vcf_cyvcf2(self, vcf_file: str) -> list[VariantEntry]:
        """
        Load variants from VCF using cyvcf2 (fast, C-based parser).

        This is 10-100x faster than pure Python parsing.
        """
        logger.info(f"Loading variants file with cyvcf2: {vcf_file}")
        variants = []

        try:
            vcf = VCF(vcf_file)

            for variant in vcf:
                chrom = variant.CHROM
                ref = variant.REF
                alt = variant.ALT[0] if variant.ALT else ""

                # Detect variant type and apply correct coordinate conversion
                is_snv = (len(ref) == len(alt) == 1)
                is_insertion = (len(ref) == 1 and len(alt) > 1)
                is_deletion = (len(ref) > 1 and len(alt) == 1)

                if is_insertion:
                    pos = variant.POS  # No conversion for insertions
                    end_pos = pos  # Insertion ends at same position
                else:
                    # SNVs and deletions: convert to 0-indexed
                    original_pos = variant.POS  # Store original 1-indexed coordinate
                    pos = original_pos - 1
                    bam_end_pos = pos + len(ref) - 1

                entry = VariantEntry(
                    chrom=chrom,
                    pos=pos,
                    end_pos=end_pos,
                    ref=ref,
                    alt=alt,
                    snp=is_snv,
                    insertion=is_insertion,
                    deletion=is_deletion,
                )
                variants.append(entry)

            vcf.close()

        except Exception as e:
            logger.error(f"Error loading VCF with cyvcf2: {e}")
            logger.info("Falling back to pure Python VCF parser")
            return self._load_vcf_python(vcf_file)

        logger.info(f"{len(variants)} variants loaded from file: {vcf_file}")
        return variants

    def _load_vcf_python(self, vcf_file: str) -> list[VariantEntry]:
        """
        Load variants from VCF using pure Python (slower but always works).

        This is the fallback method when cyvcf2 is not available.
        """
        logger.info(f"Loading variants file with Python parser: {vcf_file}")
        variants = []

        with open(vcf_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) < 5:
                    logger.error(f"Incorrectly formatted VCF entry: {line}")
                    continue

                chrom = fields[0]
                ref = fields[3]
                alt = fields[4]

                # Detect variant type and apply correct coordinate conversion
                is_snv = (len(ref) == len(alt) == 1)
                is_insertion = (len(ref) == 1 and len(alt) > 1)
                is_deletion = (len(ref) > 1 and len(alt) == 1)

                if is_insertion:
                    pos = int(fields[1])  # No conversion for insertions
                    end_pos = pos  # Insertion ends at same position
                else:
                    # SNVs and deletions: convert to 0-indexed
                    pos = int(fields[1]) - 1
                    bam_end_pos = pos + len(ref) - 1
                deletion = len(alt) == 1 and len(alt) < len(ref)

                variant = VariantEntry(
                    chrom=chrom,
                    pos=pos,
                    end_pos=bam_end_pos,
                    ref=ref,
                    alt=alt,
                    snp=is_snv,
                    insertion=is_insertion,
                    deletion=is_deletion,
                )
                variants.append(variant)

        logger.info(f"{len(variants)} variants loaded from file: {vcf_file}")
        return variants

    def load_maf(self, maf_file: str) -> list[VariantEntry]:
        """Load variants from MAF file."""
        logger.info(f"Loading variants file: {maf_file}")
        variants = []
        header_map = {}

        with open(maf_file) as f:
            # Find header line
            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    # This is the header
                    headers = line.split("\t")
                    header_map = {h: i for i, h in enumerate(headers)}
                    break

            # Validate required columns
            # Required columns for sample-agnostic workflow
            required_cols = [
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ]
            missing_cols = [col for col in required_cols if col not in header_map]
            if missing_cols:
                logger.error(f"Missing required MAF columns: {missing_cols}")
                raise ValueError("Incorrect MAF file header")

            # Load variants
            for line in f:
                line = line.strip()
                if not line:
                    continue

                fields = line.split("\t")
                variant = self._parse_maf_line(fields, header_map, line)
                if variant:
                    variants.append(variant)

        logger.info(f"{len(variants)} variants loaded from file: {maf_file}")
        return variants

    def _parse_maf_line(
        self, fields: list[str], header_map: dict[str, int], original_line: str
    ) -> VariantEntry | None:
        """Parse a single MAF line into VariantEntry."""
        try:
            gene = fields[header_map.get("Hugo_Symbol", "")]
            chrom = fields[header_map.get("Chromosome", "")]
            normalized_chrom = self._normalize_chromosome(chrom)
            original_pos = int(fields[header_map.get("Start_Position", "")])  # Store original 1-indexed
            # pos will be set by variant-type-aware logic below
            ref = fields[header_map.get("Reference_Allele", "")]
            alt = fields[header_map.get("Tumor_Seq_Allele2", "")]  # Alt should be in Allele2
            effect = fields[header_map.get("Variant_Classification", "")]
            t_ref_count = int(fields[header_map.get("t_ref_count", 0)])
            t_alt_count = int(fields[header_map.get("t_alt_count", 0)])
            # Sample-agnostic defaults
            tumor_sample = ""  # Not needed for sample-agnostic
            normal_sample = ""  # Not needed for sample-agnostic
            n_ref_count = 0  # Not needed for sample-agnostic
            n_alt_count = 0  # Not needed for sample-agnostic
            caller = ""  # Not needed for sample-agnostic

            # Detect variant type and apply correct coordinate conversion
            is_snv = (len(ref) == len(alt) == 1)
            is_insertion = (len(ref) == 1 and len(alt) > 1)
            is_deletion = (len(ref) > 1 and len(alt) == 1)

            if is_insertion:
                pos = int(fields[header_map["Start_Position"]])  # No conversion for insertions
                end_pos = pos  # Insertion ends at same position
            elif is_deletion:
                # Deletions: start at original_pos - 1, end at original_pos + len(alt) - 1
                original_pos = int(fields[header_map["Start_Position"]])
                pos = original_pos - 1  # Convert to 0-indexed
                end_pos = original_pos + len(alt) - 1  # End position for deletion
            else:
                # SNVs: convert to 0-indexed
                pos = int(fields[header_map["Start_Position"]]) - 1
            original_end_pos = int(fields[header_map["End_Position"]])  # Store original end position
            n_alt_count = int(fields[header_map["n_alt_count"]])
            effect = fields[header_map["Variant_Classification"]]

            caller = ""
            if "Caller" in header_map and len(fields) > header_map["Caller"]:
                caller = fields[header_map["Caller"]]

            # Store original MAF coordinates
            original_pos = original_pos
            original_end_pos = original_end_pos
            maf_ref = ref
            maf_alt = alt

            # Ensure end_pos is defined (for SNVs it's the same as pos)
            if not 'end_pos' in locals():
                end_pos = pos

            # Convert MAF format to VCF format
            if ref == "-":  # Insertion in MAF format
                if self.reference_getter:
                    prev_base = self.reference_getter(chrom, pos)
                    ref = prev_base
                    alt = prev_base + alt
                    end_pos -= 1
                else:
                    logger.warning(f"Cannot convert MAF insertion without reference: {chrom}:{pos}")
                    return None
            elif alt == "-":  # Deletion in MAF format
                pos -= 1
                if self.reference_getter:
                    prev_base = self.reference_getter(chrom, pos)
                    ref = prev_base + ref
                    alt = prev_base
                else:
                    logger.warning(f"Cannot convert MAF deletion without reference: {chrom}:{pos}")
                    return None
            elif len(alt) != len(ref):  # Complex indel
                pos -= 1
                if self.reference_getter:
                    prev_base = self.reference_getter(chrom, pos)
                    ref = prev_base + ref
                    alt = prev_base + alt
                else:
                    logger.warning(
                        f"Cannot convert MAF complex indel without reference: {chrom}:{pos}"
                    )
                    return None

            # Determine variant type
            snp = len(ref) == len(alt) == 1
            insertion = len(ref) == 1 and len(alt) > len(ref)
            deletion = len(alt) == 1 and len(alt) < len(ref)

            variant = VariantEntry(
                chrom=normalized_chrom,
                pos=pos,
                end_pos=end_pos,
                ref=ref,
                alt=alt,
                snp=snp,
                insertion=insertion,
                deletion=deletion,
                tumor_sample=tumor_sample,
                normal_sample=normal_sample,
                gene=gene,
                effect=effect,
                t_ref_count=t_ref_count,
                t_alt_count=t_alt_count,
                n_ref_count=n_ref_count,
                n_alt_count=n_alt_count,
                original_pos=original_pos,
                original_end_pos=original_end_pos,
                maf_ref=maf_ref,
                maf_alt=maf_alt,
                caller=caller,
                maf_line=original_line,
            )

            return variant

        except Exception as e:
            logger.error(f"Error parsing MAF line: {e}")
            return None
