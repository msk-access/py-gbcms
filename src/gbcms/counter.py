"""
Base counting algorithms for variants - Pure Python Implementation.

This module provides the standard (non-optimized) base counting implementation.
It works directly with pysam objects and is easier to debug and modify.

**When to use this module:**
- Small datasets (<10K variants)
- Development and debugging
- When Numba is not available
- When you need to modify counting logic

**Performance:** Baseline (1x)

**Alternative:** For production workloads with large datasets, see `numba_counter.py`
which provides 50-100x speedup through JIT compilation.

**Key Classes:**
- BaseCounter: Main counting class with methods for SNP, DNP, and indel variants

**Usage:**
    from gbcms.counter import BaseCounter
    from gbcms.config import Config

    config = Config(...)
    counter = BaseCounter(config)
    counter.count_bases_snp(variant, alignments, sample_name)
"""

import logging

import numpy as np
import pysam

from .config import CountType
from .variant import VariantEntry

logger = logging.getLogger(__name__)


class BaseCounter:
    """
    Performs base counting for variants using pure Python.

    This is the standard (non-optimized) implementation that works directly
    with pysam AlignedSegment objects. It's flexible and easy to debug but
    slower than the Numba-optimized version.

    **Performance Characteristics:**
    - Speed: Baseline (1x)
    - Memory: Low
    - Flexibility: High (easy to modify)
    - Debugging: Easy (pure Python)

    **For better performance on large datasets, see:**
    - `numba_counter.py` for 50-100x speedup
    - Use with multi-threading for parallel processing

    **Attributes:**
        config: Configuration object with quality thresholds and filters
        warning_counts: Track warnings to avoid spam
    """

    def __init__(self, config):
        """Initialize BaseCounter with configuration."""
        self.config = config
        self.warning_counts = {}
        # Fragment orientation tracking for Majority Rule approach
        self.fragment_orientations = {}  # frag_name -> {read1: bool, read2: bool}
        self.fragment_alleles = {}  # frag_name -> {"ref": bool, "alt": bool}
        self.processed_fragments = set()  # Track counted fragments per variant

    def _should_filter_alignment(self, aln: pysam.AlignedSegment) -> bool:
        """
        Check if alignment should be filtered based on configuration.

        Args:
            aln: BAM alignment

        Returns:
            True if alignment should be filtered (excluded)
        """
        if self.config.filter_duplicate and aln.is_duplicate:
            return True
        if self.config.filter_improper_pair and not aln.is_proper_pair:
            return True
        if self.config.filter_qc_failed and aln.is_qcfail:
            return True
        if self.config.filter_non_primary and aln.is_secondary:
            return True
        if self.config.filter_non_primary and aln.is_supplementary:
            return True
        if aln.mapping_quality < self.config.mapping_quality_threshold:
            return True
        if self.config.filter_indel and self._has_indel(aln):
            return True
        return False

    @staticmethod
    def _has_indel(aln: pysam.AlignedSegment) -> bool:
        """Check if alignment has indels."""
        if aln.cigartuples is None:
            return False
        for op, _length in aln.cigartuples:
            if op in (1, 2):  # Insertion or deletion
                return True
        return False

    def _initialize_counts(self) -> np.ndarray:
        """Initialize counts array."""
        return np.zeros(len(CountType), dtype=np.float32)

    def _initialize_fragment_maps(self) -> tuple[dict, dict, dict]:
        """Initialize fragment tracking maps."""
        return {}, {}, {}

    def _reset_fragment_tracking(self):
        """Reset fragment tracking for new variant."""
        self.fragment_orientations = {}
        self.fragment_alleles = {}
        self.processed_fragments = set()

    def _determine_fragment_orientation(self, frag_name: str) -> bool:
        """Determine fragment orientation using majority rule."""
        if frag_name not in self.fragment_orientations:
            return None  # type: ignore

        orientations = self.fragment_orientations[frag_name]

        # If we have both reads, use majority rule
        if len(orientations) == 2:
            read1_orientation = orientations.get(1, None)
            read2_orientation = orientations.get(2, None)

            if read1_orientation is not None and read2_orientation is not None:
                if read1_orientation == read2_orientation:
                    return read1_orientation  # type: ignore
                else:
                    return read1_orientation  # type: ignore

        # If we only have one read, use it
        elif len(orientations) == 1:
            return list(orientations.values())[0]  # type: ignore

        return None  # type: ignore

    def _should_count_fragment(self, frag_name: str) -> bool:
        """Check if fragment should be counted."""
        return (
            frag_name in self.fragment_orientations
            and frag_name in self.fragment_alleles
            and frag_name not in self.processed_fragments
        )

    def _track_fragment_orientation(
        self, aln: pysam.AlignedSegment, allele_type: str | None = None
    ):
        """Track fragment orientation and allele information."""
        if aln.query_name is None:
            return

        frag_name = aln.query_name
        end_no = 1 if aln.is_read1 else 2
        orientation = not aln.is_reverse

        # Initialize tracking structures
        if frag_name not in self.fragment_orientations:
            self.fragment_orientations[frag_name] = {}
        if frag_name not in self.fragment_alleles:
            self.fragment_alleles[frag_name] = {"ref": False, "alt": False}

        # Track orientation
        self.fragment_orientations[frag_name][end_no] = orientation

        # Track allele information
        if allele_type == "ref":
            self.fragment_alleles[frag_name]["ref"] = True
        elif allele_type == "alt":
            self.fragment_alleles[frag_name]["alt"] = True

    def _store_sample_counts(
        self, variant: VariantEntry, counts: np.ndarray, sample_name: str
    ) -> None:
        """Store counts for a sample, handling existing counts."""
        if sample_name not in variant.base_count:
            variant.base_count[sample_name] = counts
        else:
            variant.base_count[sample_name] += counts

    def _count_read_alleles(
        self, aln: pysam.AlignedSegment, variant: VariantEntry, counts: np.ndarray, read_pos: int
    ) -> None:
        """Count read alleles with strand analysis."""
        if aln.query_sequence is None or aln.query_qualities is None:
            return

        base = aln.query_sequence[read_pos].upper()
        qual = aln.query_qualities[read_pos]

        if qual < self.config.base_quality_threshold:
            return

        # Count total depth and strand-specific depths
        counts[CountType.DP] += 1
        if not aln.is_reverse:
            counts[CountType.DP_FORWARD] += 1
        else:
            counts[CountType.DP_REVERSE] += 1

        # Count reference/alternate alleles with strand separation
        if base == variant.ref.upper():
            counts[CountType.RD] += 1
            if not aln.is_reverse:
                counts[CountType.RD_FORWARD] += 1
            else:
                counts[CountType.RD_REVERSE] += 1
        elif base == variant.alt.upper():
            counts[CountType.AD] += 1
            if not aln.is_reverse:
                counts[CountType.AD_FORWARD] += 1
            else:
                counts[CountType.AD_REVERSE] += 1

    def _calculate_fragment_counts_basic(
        self, dpf_map: dict, rdf_map: dict, adf_map: dict
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate basic fragment counts with orientation support from tracking maps."""
        dpf = len(dpf_map)

        fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 1.0
        fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 1.0

        rdf = 0.0
        adf = 0.0

        for frag_name, end_counts in dpf_map.items():
            # Check for overlapping multimapped reads
            if any(count > 1 for count in end_counts.values()):
                if self.warning_counts["overlapping_multimap"] < self.config.max_warning_per_type:
                    logger.warning(
                        f"Fragment {frag_name} has overlapping multiple mapped alignment "
                        "at site, and will not be used"
                    )
                    self.warning_counts["overlapping_multimap"] += 1
                continue

            has_ref = frag_name in rdf_map
            has_alt = frag_name in adf_map

            if has_ref and has_alt:
                rdf += fragment_ref_weight
                adf += fragment_alt_weight
            elif has_ref:
                rdf += 1.0
            elif has_alt:
                adf += 1.0

        # Calculate orientation-specific fragment counts
        rdf_forward = sum(sum(ends.values()) for ends in rdf_map.values() if 1 in ends)
        rdf_reverse = sum(sum(ends.values()) for ends in rdf_map.values() if 2 in ends)
        adf_forward = sum(sum(ends.values()) for ends in adf_map.values() if 1 in ends)
        adf_reverse = sum(sum(ends.values()) for ends in adf_map.values() if 2 in ends)

        return dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse

    def _calculate_fragment_counts_with_orientation(
        self, dpf_map: dict, rdf_map: dict, adf_map: dict
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate fragment counts with proper orientation tracking using Majority Rule."""
        dpf = 0
        rdf = 0.0
        adf = 0.0
        rdf_forward = 0.0
        adf_forward = 0.0
        rdf_reverse = 0.0
        adf_reverse = 0.0
        fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 1.0
        fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 1.0

        # Process each fragment that has complete information
        for frag_name in list(self.fragment_orientations.keys()):
            if not self._should_count_fragment(frag_name):
                continue

            # Mark as processed
            self.processed_fragments.add(frag_name)
            dpf += 1

            # Get fragment orientation
            orientation = self._determine_fragment_orientation(frag_name)
            if orientation is None:
                continue

            # Count based on alleles present
            has_ref = self.fragment_alleles[frag_name]["ref"]
            has_alt = self.fragment_alleles[frag_name]["alt"]

            if has_ref and has_alt:
                # Fragment has both alleles
                rdf += float(fragment_ref_weight)
                adf += float(fragment_alt_weight)
                if orientation:
                    rdf_forward += float(fragment_ref_weight)
                    adf_forward += float(fragment_alt_weight)
                else:
                    rdf_reverse += float(fragment_ref_weight)
                    adf_reverse += float(fragment_alt_weight)
            elif has_ref:
                # Only reference allele
                rdf += 1.0
                if orientation:
                    rdf_forward += 1
                else:
                    rdf_reverse += 1
            elif has_alt:
                # Only alternate allele
                adf += 1.0
                if orientation:
                    adf_forward += 1
                else:
                    adf_reverse += 1

        return dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse

    def _track_fragment_indel(
        self,
        aln: pysam.AlignedSegment,
        dpf_map: dict,
        rdf_map: dict,
        adf_map: dict,
        allele_type: str,
    ) -> None:
        """Track fragment for INDEL counting with orientation support."""
        if aln.query_name is None:
            return

        end_no = 1 if aln.is_read1 else 2
        frag_name = aln.query_name

        # Initialize fragment maps with orientation tracking
        if frag_name not in dpf_map:
            dpf_map[frag_name] = {}
        if frag_name not in rdf_map:
            rdf_map[frag_name] = {}
        if frag_name not in adf_map:
            adf_map[frag_name] = {}

        # Count fragment
        dpf_map[frag_name][end_no] = dpf_map[frag_name].get(end_no, 0) + 1

        # Count allele-specific fragments with orientation
        if allele_type == "ref":
            rdf_map[frag_name][end_no] = rdf_map[frag_name].get(end_no, 0) + 1
        elif allele_type == "alt":
            adf_map[frag_name][end_no] = adf_map[frag_name].get(end_no, 0) + 1

    def _calculate_fragment_counts_indel(
        self, dpf_map: dict, rdf_map: dict, adf_map: dict
    ) -> tuple[float, float, float, float, float, float, float]:
        """Calculate fragment counts for INDEL from tracking maps with orientation support."""
        dpf = len(dpf_map)

        fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 1.0
        fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 1.0

        rdf = 0.0
        adf = 0.0

        for frag_name, end_counts in dpf_map.items():
            # Check for overlapping multimapped reads
            if any(count > 1 for count in end_counts.values()):
                if self.warning_counts["overlapping_multimap"] < self.config.max_warning_per_type:
                    logger.warning(
                        f"Fragment {frag_name} has overlapping multiple mapped alignment "
                        "at site, and will not be used"
                    )
                    self.warning_counts["overlapping_multimap"] += 1
                continue

            has_ref = frag_name in rdf_map
            has_alt = frag_name in adf_map

            if has_ref and has_alt:
                rdf += fragment_ref_weight
                adf += fragment_alt_weight
            elif has_ref:
                rdf += 1.0
            elif has_alt:
                adf += 1.0

        # Calculate orientation-specific fragment counts
        rdf_forward = sum(sum(ends.values()) for ends in rdf_map.values() if 1 in ends)
        rdf_reverse = sum(sum(ends.values()) for ends in rdf_map.values() if 2 in ends)
        adf_forward = sum(sum(ends.values()) for ends in adf_map.values() if 1 in ends)
        adf_reverse = sum(sum(ends.values()) for ends in adf_map.values() if 2 in ends)

        return dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse

    def _track_fragment_generic(
        self,
        aln: pysam.AlignedSegment,
        dpf_map: dict,
        rdf_map: dict,
        adf_map: dict,
        base: str,
        variant: VariantEntry,
    ) -> None:
        """Track fragment for generic counting with orientation support."""
        if aln.query_name is None:
            return

        end_no = 1 if aln.is_read1 else 2
        frag_name = aln.query_name

        # Initialize fragment maps with orientation tracking
        if frag_name not in dpf_map:
            dpf_map[frag_name] = {}
        if frag_name not in rdf_map:
            rdf_map[frag_name] = {}
        if frag_name not in adf_map:
            adf_map[frag_name] = {}

        # Count fragment
        dpf_map[frag_name][end_no] = dpf_map[frag_name].get(end_no, 0) + 1

        # Count allele-specific fragments with orientation
        if base == variant.ref.upper():
            rdf_map[frag_name][end_no] = rdf_map[frag_name].get(end_no, 0) + 1
        elif base == variant.alt.upper():
            adf_map[frag_name][end_no] = adf_map[frag_name].get(end_no, 0) + 1

    def count_bases_snp(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for SNP variants.

        This method processes alignments overlapping a SNP variant position and
        counts reference and alternate alleles, considering strand orientation
        and fragment relationships when enabled.

        Args:
            variant: VariantEntry object containing position and alleles
            alignments: List of pysam AlignedSegment objects overlapping variant
            sample_name: Name of sample for storing results

        Note:
            Applies all configured filters before counting
            Updates variant.base_count[sample_name] with results
        """
        counts = self._initialize_counts()

        # Fragment tracking for fragment counts
        dpf_map, rdf_map, adf_map = self._initialize_fragment_maps()

        # Reset fragment tracking for new variant
        self._reset_fragment_tracking()

        for aln in alignments:
            if self._should_filter_alignment(aln):
                continue

            # Check if alignment overlaps variant position
            if (aln.reference_start is not None and aln.reference_start > variant.pos) or (
                aln.reference_end is not None and aln.reference_end <= variant.pos
            ):
                continue

            # Get the base at variant position
            read_pos = None
            for read_idx, ref_idx in aln.get_aligned_pairs(matches_only=False):
                if ref_idx == variant.pos:
                    read_pos = read_idx
                    break

            if read_pos is None:
                continue  # Variant position is in deletion

            # Use helper method for allele counting
            self._count_read_alleles(aln, variant, counts, read_pos)

            # Track fragments if enabled
            if self.config.output_fragment_count and aln.query_name is not None:
                # Track orientation and allele information for this fragment
                base = aln.query_sequence[read_pos].upper() if aln.query_sequence else None
                allele_type = None
                if base == variant.ref.upper():
                    allele_type = "ref"
                elif base == variant.alt.upper():
                    allele_type = "alt"

                self._track_fragment_orientation(aln, allele_type)

        # Calculate fragment counts using helper
        if self.config.output_fragment_count:
            dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse = (
                self._calculate_fragment_counts_with_orientation(dpf_map, rdf_map, adf_map)
            )
            counts[CountType.DPF] = dpf
            counts[CountType.RDF] = rdf
            counts[CountType.ADF] = adf
            counts[CountType.RDF_FORWARD] = rdf_forward
            counts[CountType.RDF_REVERSE] = rdf_reverse
            counts[CountType.ADF_FORWARD] = adf_forward
            counts[CountType.ADF_REVERSE] = adf_reverse

        # Store counts using helper
        self._store_sample_counts(variant, counts, sample_name)

    def count_bases_indel(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for insertion/deletion variants.

        Handles complex indel variants by examining CIGAR strings and
        sequence alignments to accurately count reference and alternate alleles.

        Args:
            variant: VariantEntry with indel information
            alignments: List of overlapping alignments
            sample_name: Sample name for results storage
        """
        counts = self._initialize_counts()

        # Fragment tracking for fragment counts
        dpf_map, rdf_map, adf_map = self._initialize_fragment_maps()

        # Reset fragment tracking for new variant
        self._reset_fragment_tracking()

        for aln in alignments:
            if self._should_filter_alignment(aln):
                continue

            # Parse CIGAR to find indels at the variant position
            ref_pos = aln.reference_start
            read_pos = 0
            if aln.cigartuples is None:
                continue

            for i, (cigar_op, cigar_len) in enumerate(aln.cigartuples):
                if ref_pos is not None and ref_pos > variant.pos + 1:
                    break

                if cigar_op == 0:  # Match/mismatch (M)
                    if ref_pos is not None and ref_pos + cigar_len - 1 == variant.pos:
                        # Found the variant position in a match
                        if i + 1 < len(aln.cigartuples):
                            next_op, next_len = aln.cigartuples[i + 1]
                            if next_op == 1 and variant.insertion:  # Insertion (I)
                                # Handle insertion case
                                read_pos += cigar_len
                                # Check for base at insertion position
                                if (
                                    aln.query_sequence is not None
                                    and read_pos < len(aln.query_sequence)
                                    and aln.query_sequence[read_pos].upper() == variant.alt.upper()
                                ):
                                    # Use helper for allele counting
                                    self._count_read_alleles(aln, variant, counts, read_pos)
                                    # Track fragments if enabled
                                    if (
                                        self.config.output_fragment_count
                                        and aln.query_name is not None
                                    ):
                                        self._track_fragment_orientation(aln, "alt")
                                    break
                            elif next_op == 2 and variant.deletion:  # Deletion (D)
                                # Handle deletion case
                                if (
                                    ref_pos is not None
                                    and ref_pos <= variant.pos + 1 < ref_pos + cigar_len
                                ):
                                    offset = variant.pos + 1 - ref_pos
                                    # For deletions, we count the reference base
                                    if (
                                        aln.query_sequence is not None
                                        and read_pos + offset < len(aln.query_sequence)
                                        and aln.query_sequence[read_pos + offset].upper()
                                        == variant.ref.upper()
                                    ):
                                        self._count_read_alleles(
                                            aln, variant, counts, read_pos + offset
                                        )
                                        # Track fragments if enabled
                                        if (
                                            self.config.output_fragment_count
                                            and aln.query_name is not None
                                        ):
                                            self._track_fragment_orientation(aln, "ref")
                                        break

                    if ref_pos is not None:
                        ref_pos += cigar_len
                    read_pos += cigar_len

                elif cigar_op == 1:  # Insertion (I)
                    read_pos += cigar_len
                elif cigar_op == 2:  # Deletion (D)
                    if ref_pos is not None:
                        ref_pos += cigar_len
                elif cigar_op == 4:  # Soft clip (S)
                    read_pos += cigar_len
                elif cigar_op == 5:  # Hard clip (H)
                    pass
                elif cigar_op == 3:  # Skipped region (N)
                    if ref_pos is not None:
                        ref_pos += cigar_len

        # Calculate fragment counts inline (no separate traversal)
        if self.config.output_fragment_count:
            dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse = (
                self._calculate_fragment_counts_with_orientation(dpf_map, rdf_map, adf_map)
            )
            counts[CountType.DPF] = dpf
            counts[CountType.RDF] = rdf
            counts[CountType.ADF] = adf
            counts[CountType.RDF_FORWARD] = rdf_forward
            counts[CountType.RDF_REVERSE] = rdf_reverse
            counts[CountType.ADF_FORWARD] = adf_forward
            counts[CountType.ADF_REVERSE] = adf_reverse

        # Store counts using helper
        self._store_sample_counts(variant, counts, sample_name)

    def count_bases_generic(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Complete generic counting with full feature support for all variant types.

        Provides comprehensive analysis including:
        - Read counts (DP, RD, AD)
        - Complete strand analysis (forward/reverse)
        - Fragment counting (DPF, RDF, ADF)
        - Statistical analysis support (strand bias)
        - Complete filtering (all 7 conditions)
        - CIGAR string parsing for complex variants

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts

        Raises:
            ValueError: If inputs are invalid or counting fails
        """
        try:
            # Validate inputs
            if variant is None:
                raise ValueError("Variant cannot be None")
            if not alignments:
                raise ValueError(
                    f"No alignments provided for variant at {variant.get_variant_key()}"
                )
            if not sample_name:
                raise ValueError("Sample name cannot be empty")

            counts = self._initialize_counts()

            # Fragment tracking for fragment counts
            dpf_map, rdf_map, adf_map = self._initialize_fragment_maps()

            # Reset fragment tracking for new variant
            self._reset_fragment_tracking()

            # Process each alignment with complete filtering and counting
            for aln in alignments:
                # Apply all 7 filter conditions
                if self._should_filter_alignment(aln):
                    continue

                # Find the base at variant position using aligned pairs
                base = None
                read_idx = None

                for read_idx, ref_idx in aln.get_aligned_pairs(matches_only=False):
                    if ref_idx == variant.pos and read_idx is not None:
                        if (
                            aln.query_sequence is not None
                            and aln.query_qualities is not None
                            and read_idx < len(aln.query_sequence)
                            and read_idx < len(aln.query_qualities)
                        ):
                            # Check base quality threshold
                            if aln.query_qualities[read_idx] >= self.config.base_quality_threshold:
                                base = aln.query_sequence[read_idx].upper()
                        break

                if base is None or read_idx is None:
                    continue

                # Use helper method for allele counting
                self._count_read_alleles(aln, variant, counts, read_idx)

                # Track fragments if enabled (single traversal)
                if self.config.output_fragment_count and aln.query_name is not None:
                    # Track orientation and allele information for this fragment
                    allele_type = None
                    if base == variant.ref.upper():
                        allele_type = "ref"
                    elif base == variant.alt.upper():
                        allele_type = "alt"

                    self._track_fragment_orientation(aln, allele_type)

            # Calculate fragment counts inline (no separate traversal)
            if self.config.output_fragment_count:
                dpf, rdf, adf, rdf_forward, rdf_reverse, adf_forward, adf_reverse = (
                    self._calculate_fragment_counts_with_orientation(dpf_map, rdf_map, adf_map)
                )
                counts[CountType.DPF] = dpf
                counts[CountType.RDF] = rdf
                counts[CountType.ADF] = adf
                counts[CountType.RDF_FORWARD] = rdf_forward
                counts[CountType.RDF_REVERSE] = rdf_reverse
                counts[CountType.ADF_FORWARD] = adf_forward
                counts[CountType.ADF_REVERSE] = adf_reverse

            # Store counts using helper
            self._store_sample_counts(variant, counts, sample_name)

        except AttributeError as e:
            raise ValueError(f"Invalid variant object: {e}") from None
        except Exception as e:
            variant_key = variant.get_variant_key() if variant else "unknown"
            raise ValueError(f"Counting failed for variant {variant_key}: {e}") from e

    def calculate_strand_bias(
        self,
        ref_forward: int,
        ref_reverse: int,
        alt_forward: int,
        alt_reverse: int,
        min_depth: int = 10,
    ) -> tuple[float, float, str]:
        """
        Calculate strand bias using Fisher's exact test.

        This unified method works with any orientation-specific counts:
        - Read-level: RD_FORWARD, RD_REVERSE, AD_FORWARD, AD_REVERSE
        - Fragment-level: RDF_FORWARD, RDF_REVERSE, ADF_FORWARD, ADF_REVERSE
        """
        try:
            import numpy as np
            from scipy.stats import fisher_exact

            # Check minimum depth requirement
            total_depth = ref_forward + ref_reverse + alt_forward + alt_reverse
            if total_depth < min_depth:
                return 1.0, 1.0, "insufficient_depth"

            # Create 2x2 contingency table
            # [[ref_forward, ref_reverse],
            #  [alt_forward, alt_reverse]]
            table = np.array([[ref_forward, ref_reverse], [alt_forward, alt_reverse]])

            # Fisher's exact test
            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

            # Determine bias direction
            total_forward = ref_forward + alt_forward
            total_reverse = ref_reverse + alt_reverse

            if total_forward > 0 and total_reverse > 0:
                forward_ratio = ref_forward / total_forward if total_forward > 0 else 0
                reverse_ratio = ref_reverse / total_reverse if total_reverse > 0 else 0

                if forward_ratio > reverse_ratio + 0.1:  # 10% threshold
                    bias_direction = "forward"
                elif reverse_ratio > forward_ratio + 0.1:
                    bias_direction = "reverse"
                else:
                    bias_direction = "none"
            else:
                bias_direction = "none"

            return p_value, odds_ratio, bias_direction

        except ImportError:
            logger.warning("scipy not available for strand bias calculation")
            return 1.0, 1.0, "scipy_unavailable"
        except Exception as e:
            logger.warning(f"Error calculating strand bias: {e}")
            return 1.0, 1.0, "error"

    def smart_count_variant(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count variant bases using standard counting for SNV/Insertion/Deletion, generic for Complex.

        Routing strategy:
        - SNP variants → count_bases_snp() (standard SNP counting)
        - Insertion variants → count_bases_indel() (standard indel counting)
        - Deletion variants → count_bases_indel() (standard indel counting)
        - Complex variants → count_bases_generic() (generic counting with full features)

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts

        Raises:
            ValueError: If variant type cannot be determined or counting fails
        """
        try:
            # Validate inputs
            if variant is None:
                raise ValueError("Variant cannot be None")
            if not alignments:
                raise ValueError(
                    f"No alignments provided for variant at {variant.get_variant_key()}"
                )
            if not sample_name:
                raise ValueError("Sample name cannot be empty")

            # Use standard counting for SNV, Insertion, Deletion; generic for Complex
            if variant.insertion or variant.deletion:
                self.count_bases_indel(variant, alignments, sample_name)
            elif variant.snp:
                self.count_bases_snp(variant, alignments, sample_name)
            else:
                # Complex variants (including previous DNP) use generic counting
                self.count_bases_generic(variant, alignments, sample_name)

        except AttributeError as e:
            raise ValueError(f"Invalid variant object: {e}") from None
        except Exception as e:
            variant_key = variant.get_variant_key() if variant else "unknown"
            raise ValueError(f"Counting failed for variant {variant_key}: {e}") from e
