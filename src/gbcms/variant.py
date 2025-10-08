"""Variant loading and representation."""

import logging
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


@dataclass
class VariantEntry:
    """Represents a variant with its counts across samples."""

    chrom: str
    pos: int  # 0-indexed
    end_pos: int
    ref: str
    alt: str
    snp: bool = False
    dnp: bool = False
    dnp_len: int = 0
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
    maf_pos: int = 0
    maf_end_pos: int = 0
    maf_ref: str = ""
    maf_alt: str = ""
    caller: str = ""
    base_count: dict[str, np.ndarray] = field(default_factory=dict)
    duplicate_variant_ptr: Optional["VariantEntry"] = None
    maf_line: str = ""  # Store original MAF line for output

    def get_variant_key(self) -> tuple[str, int, str, str]:
        """Return unique key for variant identification."""
        return (self.chrom, self.pos, self.ref, self.alt)

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
        return self.pos < other.pos

    def _chrom_sort_key(self) -> tuple:
        """Generate sort key for chromosome."""
        chrom = self.chrom.replace("chr", "")
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

    def __init__(self, reference_getter=None):
        """
        Initialize variant loader.

        Args:
            reference_getter: Callable that takes (chrom, pos) and returns base
        """
        self.reference_getter = reference_getter

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
                pos = variant.POS - 1  # Convert to 0-indexed
                ref = variant.REF

                # Handle multiple alts - take first one
                alt = variant.ALT[0] if variant.ALT else ""

                end_pos = pos + len(ref) - 1

                # Determine variant type
                snp = len(alt) == len(ref) == 1
                dnp = len(alt) == len(ref) > 1
                dnp_len = len(ref) if dnp else 0
                insertion = len(alt) > len(ref)
                deletion = len(alt) < len(ref)

                entry = VariantEntry(
                    chrom=chrom,
                    pos=pos,
                    end_pos=end_pos,
                    ref=ref,
                    alt=alt,
                    snp=snp,
                    dnp=dnp,
                    dnp_len=dnp_len,
                    insertion=insertion,
                    deletion=deletion,
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
                pos = int(fields[1]) - 1  # Convert to 0-indexed
                ref = fields[3]
                alt = fields[4]

                # Handle multiple alts - take first one
                if "," in alt:
                    alt = alt.split(",")[0]

                end_pos = pos + len(ref) - 1

                # Determine variant type
                snp = len(alt) == len(ref) == 1
                dnp = len(alt) == len(ref) > 1
                dnp_len = len(ref) if dnp else 0
                insertion = len(alt) > len(ref)
                deletion = len(alt) < len(ref)

                variant = VariantEntry(
                    chrom=chrom,
                    pos=pos,
                    end_pos=end_pos,
                    ref=ref,
                    alt=alt,
                    snp=snp,
                    dnp=dnp,
                    dnp_len=dnp_len,
                    insertion=insertion,
                    deletion=deletion,
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
            required_cols = [
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele1",
                "Tumor_Seq_Allele2",
                "Tumor_Sample_Barcode",
                "Matched_Norm_Sample_Barcode",
                "t_ref_count",
                "t_alt_count",
                "n_ref_count",
                "n_alt_count",
                "Variant_Classification",
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
            gene = fields[header_map["Hugo_Symbol"]]
            chrom = fields[header_map["Chromosome"]]
            pos = int(fields[header_map["Start_Position"]]) - 1  # Convert to 0-indexed
            end_pos = int(fields[header_map["End_Position"]]) - 1
            ref = fields[header_map["Reference_Allele"]]
            alt = fields[header_map["Tumor_Seq_Allele1"]]

            # Use Tumor_Seq_Allele2 if Allele1 is empty or same as ref
            if not alt or alt == ref:
                alt = fields[header_map["Tumor_Seq_Allele2"]]

            if not alt or alt == ref:
                logger.warning(f"Could not find valid alt allele for variant: {chrom}:{pos + 1}")
                return None

            tumor_sample = fields[header_map["Tumor_Sample_Barcode"]]
            normal_sample = fields[header_map["Matched_Norm_Sample_Barcode"]]
            t_ref_count = int(fields[header_map["t_ref_count"]])
            t_alt_count = int(fields[header_map["t_alt_count"]])
            n_ref_count = int(fields[header_map["n_ref_count"]])
            n_alt_count = int(fields[header_map["n_alt_count"]])
            effect = fields[header_map["Variant_Classification"]]

            caller = ""
            if "Caller" in header_map and len(fields) > header_map["Caller"]:
                caller = fields[header_map["Caller"]]

            # Store original MAF coordinates
            maf_pos = pos
            maf_end_pos = end_pos
            maf_ref = ref
            maf_alt = alt

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
            snp = len(alt) == len(ref) == 1
            dnp = len(alt) == len(ref) > 1
            dnp_len = len(ref) if dnp else 0
            insertion = len(alt) > len(ref)
            deletion = len(alt) < len(ref)

            variant = VariantEntry(
                chrom=chrom,
                pos=pos,
                end_pos=end_pos,
                ref=ref,
                alt=alt,
                snp=snp,
                dnp=dnp,
                dnp_len=dnp_len,
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
                maf_pos=maf_pos,
                maf_end_pos=maf_end_pos,
                maf_ref=maf_ref,
                maf_alt=maf_alt,
                caller=caller,
                maf_line=original_line,
            )

            return variant

        except Exception as e:
            logger.error(f"Error parsing MAF line: {e}")
            return None
