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
from collections import defaultdict

import numpy as np
import pysam

from .config import Config, CountType
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
    - Use with `parallel.py` for multi-core processing

    **Attributes:**
        config: Configuration object with quality thresholds and filters
        warning_counts: Track warnings to avoid spam
    """

    def __init__(self, config: Config):
        """
        Initialize base counter.

        Args:
            config: Configuration object with quality filters and thresholds
        """
        self.config = config
        self.warning_counts: dict[str, int] = defaultdict(int)

    def filter_alignment(self, aln: pysam.AlignedSegment) -> bool:
        """
        Check if alignment should be filtered.

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

    def count_bases_snp(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for SNP variants.

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts
        """
        counts = np.zeros(len(CountType), dtype=np.float32)

        # Fragment tracking for fragment counts
        dpf_map: dict[str, dict[int, int]] = {}
        rdf_map: dict[str, dict[int, int]] = {}
        adf_map: dict[str, dict[int, int]] = {}

        for aln in alignments:
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

            # Check if query sequence and qualities are available
            if aln.query_sequence is None or aln.query_qualities is None:
                continue

            # Get base and quality
            base = aln.query_sequence[read_pos].upper()
            qual = aln.query_qualities[read_pos]

            if qual < self.config.base_quality_threshold:
                continue

            # Count total depth
            counts[CountType.DP] += 1
            if not aln.is_reverse:
                counts[CountType.DPP] += 1

            # Track fragment
            end_no = 1 if aln.is_read1 else 2
            if self.config.output_fragment_count:
                if aln.query_name is not None:
                    if aln.query_name not in dpf_map:
                        dpf_map[aln.query_name] = {}
                    dpf_map[aln.query_name][end_no] = dpf_map[aln.query_name].get(end_no, 0) + 1

            # Count ref/alt
            if base == variant.ref:
                counts[CountType.RD] += 1
                if not aln.is_reverse:
                    counts[CountType.RDP] += 1
                if self.config.output_fragment_count and aln.query_name is not None:
                    if aln.query_name not in rdf_map:
                        rdf_map[aln.query_name] = {}
                    rdf_map[aln.query_name][end_no] = rdf_map[aln.query_name].get(end_no, 0) + 1
            elif base == variant.alt:
                counts[CountType.AD] += 1
                if not aln.is_reverse:
                    counts[CountType.ADP] += 1
                if self.config.output_fragment_count and aln.query_name is not None:
                    if aln.query_name not in adf_map:
                        adf_map[aln.query_name] = {}
                    adf_map[aln.query_name][end_no] = adf_map[aln.query_name].get(end_no, 0) + 1

        # Calculate strand bias for this sample
        ref_forward = int(counts[CountType.RDP])
        ref_reverse = int(counts[CountType.RD]) - ref_forward
        alt_forward = int(counts[CountType.ADP])
        alt_reverse = int(counts[CountType.AD]) - alt_forward

        strand_bias_pval, strand_bias_or, strand_bias_dir = self.calculate_strand_bias(
            ref_forward, ref_reverse, alt_forward, alt_reverse
        )

        # Note: Strand bias is calculated on-the-fly during output, not stored here

        # Calculate fragment counts
        if self.config.output_fragment_count:
            counts[CountType.DPF] = len(dpf_map)

            fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 0
            fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 0

            for frag_name, end_counts in dpf_map.items():
                # Check for overlapping multimapped reads
                if any(count > 1 for count in end_counts.values()):
                    if (
                        self.warning_counts["overlapping_multimap"]
                        < self.config.max_warning_per_type
                    ):
                        logger.warning(
                            f"Fragment {frag_name} has overlapping multiple mapped alignment "
                            f"at site: {variant.chrom}:{variant.pos + 1}, and will not be used"
                        )
                        self.warning_counts["overlapping_multimap"] += 1
                    continue

                has_ref = frag_name in rdf_map
                has_alt = frag_name in adf_map

                if has_ref and has_alt:
                    counts[CountType.RDF] += fragment_ref_weight
                    counts[CountType.ADF] += fragment_alt_weight
                elif has_ref:
                    counts[CountType.RDF] += 1
                elif has_alt:
                    counts[CountType.ADF] += 1

        # Store counts
        if sample_name not in variant.base_count:
            variant.base_count[sample_name] = counts
        else:
            variant.base_count[sample_name] += counts

    def count_bases_dnp(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for DNP (di-nucleotide polymorphism) variants.

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts
        """
        counts = np.zeros(len(CountType), dtype=np.float32)

        dpf_map: dict[str, dict[int, int]] = {}
        rdf_map: dict[str, dict[int, int]] = {}
        adf_map: dict[str, dict[int, int]] = {}

        for aln in alignments:
            # Check if alignment fully covers the DNP
            if (aln.reference_start is not None and aln.reference_start > variant.pos) or (
                aln.reference_end is not None
                and aln.reference_end <= variant.pos + variant.dnp_len - 1
            ):
                continue

            # Find the read positions corresponding to the DNP
            read_bases = []
            for read_idx, ref_idx in aln.get_aligned_pairs(matches_only=True):
                if ref_idx is not None and variant.pos <= ref_idx < variant.pos + variant.dnp_len:
                    if aln.query_sequence is not None:
                        read_bases.append((read_idx, aln.query_sequence[read_idx]))

            if len(read_bases) != variant.dnp_len:
                continue  # DNP not fully covered

            # Check if query sequence and qualities are available
            if aln.query_sequence is None or aln.query_qualities is None:
                continue

            # Get the DNP sequence and minimum quality
            dnp_seq = "".join([base for _, base in read_bases]).upper()
            min_qual = min([aln.query_qualities[idx] for idx, _ in read_bases])

            if min_qual < self.config.base_quality_threshold:
                continue

            # Count total depth
            counts[CountType.DP] += 1
            if not aln.is_reverse:
                counts[CountType.DPP] += 1

            # Track fragment
            end_no = 1 if aln.is_read1 else 2
            if self.config.output_fragment_count:
                if aln.query_name is not None:
                    if aln.query_name not in dpf_map:
                        dpf_map[aln.query_name] = {}
                    dpf_map[aln.query_name][end_no] = dpf_map[aln.query_name].get(end_no, 0) + 1

            # Count ref/alt
            if dnp_seq == variant.ref:
                counts[CountType.RD] += 1
                if not aln.is_reverse:
                    counts[CountType.RDP] += 1
                if self.config.output_fragment_count and aln.query_name is not None:
                    if aln.query_name not in rdf_map:
                        rdf_map[aln.query_name] = {}
                    rdf_map[aln.query_name][end_no] = rdf_map[aln.query_name].get(end_no, 0) + 1
            elif dnp_seq == variant.alt:
                counts[CountType.AD] += 1
                if not aln.is_reverse:
                    counts[CountType.ADP] += 1
                if self.config.output_fragment_count and aln.query_name is not None:
                    if aln.query_name not in adf_map:
                        adf_map[aln.query_name] = {}
                    adf_map[aln.query_name][end_no] = adf_map[aln.query_name].get(end_no, 0) + 1

        # Calculate fragment counts
        if self.config.output_fragment_count:
            counts[CountType.DPF] = len(dpf_map)

            fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 0
            fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 0

            for frag_name, end_counts in dpf_map.items():
                if any(count > 1 for count in end_counts.values()):
                    continue

                has_ref = frag_name in rdf_map
                has_alt = frag_name in adf_map

                if has_ref and has_alt:
                    counts[CountType.RDF] += fragment_ref_weight
                    counts[CountType.ADF] += fragment_alt_weight
                elif has_ref:
                    counts[CountType.RDF] += 1
                elif has_alt:
                    counts[CountType.ADF] += 1

        if sample_name not in variant.base_count:
            variant.base_count[sample_name] = counts
        else:
            variant.base_count[sample_name] += counts

    def count_bases_indel(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for indel variants using DMP method.

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts
        """
        counts = np.zeros(len(CountType), dtype=np.float32)

        dpf_map: dict[str, dict[int, int]] = {}
        rdf_map: dict[str, dict[int, int]] = {}
        adf_map: dict[str, dict[int, int]] = {}

        for aln in alignments:
            # Check if alignment overlaps the indel region
            if (aln.reference_start is not None and aln.reference_start > variant.pos + 1) or (
                aln.reference_end is not None and aln.reference_end <= variant.pos
            ):
                continue

            # Parse CIGAR to find indels at the variant position
            matched_indel = False
            ref_pos = aln.reference_start
            read_pos = 0

            if aln.cigartuples is None:
                continue

            for i, (cigar_op, cigar_len) in enumerate(aln.cigartuples):
                if ref_pos is not None and ref_pos > variant.pos + 1:
                    break

                if cigar_op == 0:  # Match/mismatch (M)
                    # Check if variant position is at the end of this match
                    if ref_pos is not None and ref_pos + cigar_len - 1 == variant.pos:
                        # Look ahead for insertion or deletion
                        if i + 1 < len(aln.cigartuples):
                            next_op, next_len = aln.cigartuples[i + 1]

                            if next_op == 1 and variant.insertion:  # Insertion (I)
                                expected_ins_len = len(variant.alt) - len(variant.ref)
                                if next_len == expected_ins_len:
                                    # Check if insertion sequence matches
                                    if aln.query_sequence is not None:
                                        ins_seq = aln.query_sequence[
                                            read_pos + cigar_len : read_pos + cigar_len + next_len
                                        ]
                                        expected_ins_seq = variant.alt[len(variant.ref) :]
                                        if ins_seq == expected_ins_seq:
                                            matched_indel = True
                            elif next_op == 2 and variant.deletion:  # Deletion (D)
                                expected_del_len = len(variant.ref) - len(variant.alt)
                                if next_len == expected_del_len:
                                    matched_indel = True

                    # Check if we can count depth at pos+1
                    if ref_pos is not None and ref_pos <= variant.pos + 1 < ref_pos + cigar_len:
                        # Check if query qualities are available
                        if aln.query_qualities is None:
                            continue

                        # Get base quality at pos+1
                        offset = variant.pos + 1 - ref_pos
                        qual = aln.query_qualities[read_pos + offset]

                        if qual >= self.config.base_quality_threshold:
                            # Count total depth
                            counts[CountType.DP] += 1
                            if not aln.is_reverse:
                                counts[CountType.DPP] += 1

                            # Track fragment
                            end_no = 1 if aln.is_read1 else 2
                            if self.config.output_fragment_count:
                                if aln.query_name is not None:
                                    if aln.query_name not in dpf_map:
                                        dpf_map[aln.query_name] = {}
                                    dpf_map[aln.query_name][end_no] = (
                                        dpf_map[aln.query_name].get(end_no, 0) + 1
                                    )

                            # Count ref/alt based on matched indel
                            if matched_indel:
                                counts[CountType.AD] += 1
                                if not aln.is_reverse:
                                    counts[CountType.ADP] += 1
                                if self.config.output_fragment_count and aln.query_name is not None:
                                    if aln.query_name not in adf_map:
                                        adf_map[aln.query_name] = {}
                                    adf_map[aln.query_name][end_no] = (
                                        adf_map[aln.query_name].get(end_no, 0) + 1
                                    )
                            else:
                                counts[CountType.RD] += 1
                                if not aln.is_reverse:
                                    counts[CountType.RDP] += 1
                                if self.config.output_fragment_count and aln.query_name is not None:
                                    if aln.query_name not in rdf_map:
                                        rdf_map[aln.query_name] = {}
                                    rdf_map[aln.query_name][end_no] = (
                                        rdf_map[aln.query_name].get(end_no, 0) + 1
                                    )

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

        # Calculate fragment counts
        if self.config.output_fragment_count:
            counts[CountType.DPF] = len(dpf_map)

            fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 0
            fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 0

            for frag_name, end_counts in dpf_map.items():
                if any(count > 1 for count in end_counts.values()):
                    continue

                has_ref = frag_name in rdf_map
                has_alt = frag_name in adf_map

                if has_ref and has_alt:
                    counts[CountType.RDF] += fragment_ref_weight
                    counts[CountType.ADF] += fragment_alt_weight
                elif has_ref:
                    counts[CountType.RDF] += 1
                elif has_alt:
                    counts[CountType.ADF] += 1

        # Store counts
        if sample_name not in variant.base_count:
            variant.base_count[sample_name] = counts
        else:
            variant.base_count[sample_name] += counts

    def count_bases_generic(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Generic counting algorithm that works for all variant types.

        This algorithm extracts the alignment allele by parsing CIGAR and comparing
        directly to ref/alt. Works better for complex variants but may give slightly
        different results than the specialized counting methods.

        This is equivalent to the C++ baseCountGENERIC function.

        Args:
            variant: Variant to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name
        """
        counts = np.zeros(len(CountType), dtype=np.float32)
        dpf_map: dict[str, dict[int, int]] = {}
        rdf_map: dict[str, dict[int, int]] = {}
        adf_map: dict[str, dict[int, int]] = {}

        for aln in alignments:
            if self.filter_alignment(aln):
                continue

            # Check if alignment overlaps variant region
            if aln.reference_end is not None and (
                aln.reference_end <= variant.pos or aln.reference_start > variant.end_pos
            ):
                continue

            # Extract alignment allele by parsing CIGAR
            alignment_allele = ""
            cur_bq = float("inf")
            partially_cover = False

            if aln.reference_start > variant.pos or (
                aln.reference_end is not None and aln.reference_end < variant.end_pos
            ):
                partially_cover = True

            # Check if query sequence and qualities are available
            if aln.query_sequence is None or aln.query_qualities is None:
                continue

            # Parse CIGAR to extract allele
            ref_pos = aln.reference_start
            read_pos = 0
            additional_insertion = False

            if aln.cigartuples is None:
                continue

            for i, (op, length) in enumerate(aln.cigartuples):
                if (
                    aln.reference_end is not None
                    and ref_pos > variant.end_pos
                    and not additional_insertion
                ):
                    break

                if op == 0:  # M (match/mismatch)
                    if ref_pos is not None and ref_pos + length - 1 >= variant.pos:
                        start_idx = read_pos + max(variant.pos, ref_pos) - ref_pos
                        str_len = min(
                            length,
                            min(variant.end_pos, ref_pos + length - 1)
                            + 1
                            - max(variant.pos, ref_pos),
                        )
                        alignment_allele += aln.query_sequence[start_idx : start_idx + str_len]

                        # Get minimum base quality
                        for bq_idx in range(str_len):
                            cur_bq = min(cur_bq, aln.query_qualities[start_idx + bq_idx])

                    if ref_pos is not None:
                        ref_pos += length
                    read_pos += length

                    # Allow additional insertion if M falls at variant end
                    if ref_pos is not None and ref_pos == variant.end_pos + 1:
                        if i + 1 < len(aln.cigartuples) and aln.cigartuples[i + 1][0] == 1:
                            additional_insertion = True

                elif op == 1:  # I (insertion)
                    if ref_pos is not None and ref_pos >= variant.pos:
                        alignment_allele += aln.query_sequence[read_pos : read_pos + length]
                        for bq_idx in range(length):
                            cur_bq = min(cur_bq, aln.query_qualities[read_pos + bq_idx])
                    read_pos += length
                    additional_insertion = False

                elif op == 4:  # S (soft clip)
                    read_pos += length

                elif op in [2, 3]:  # D or N (deletion/skip)
                    if (
                        aln.reference_end is not None
                        and ref_pos is not None
                        and ref_pos + length - 1 > variant.end_pos
                    ):
                        alignment_allele = "U"  # Unmatched deletion
                    if ref_pos is not None:
                        ref_pos += length

                    # Allow additional insertion if D/N falls at variant end
                    if ref_pos is not None and ref_pos == variant.end_pos + 1:
                        if i + 1 < len(aln.cigartuples) and aln.cigartuples[i + 1][0] == 1:
                            additional_insertion = True

            # Check base quality threshold
            if cur_bq < self.config.base_quality_threshold:
                continue

            # Count depth
            counts[CountType.DP] += 1
            if not aln.is_reverse:
                counts[CountType.DPP] += 1

            # Track fragment
            end_no = 1 if aln.is_read1 else 2
            frag_name = aln.query_name

            if self.config.output_fragment_count:
                if frag_name is not None:
                    if frag_name not in dpf_map:
                        dpf_map[frag_name] = {}
                    if end_no not in dpf_map[frag_name]:
                        dpf_map[frag_name][end_no] = 0
                    dpf_map[frag_name][end_no] += 1

            # Count ref/alt (skip if partially covered)
            if not partially_cover:
                if alignment_allele == variant.ref:
                    counts[CountType.RD] += 1
                    if not aln.is_reverse:
                        counts[CountType.RDP] += 1

                    if self.config.output_fragment_count and frag_name is not None:
                        if frag_name not in rdf_map:
                            rdf_map[frag_name] = {}
                        if end_no not in rdf_map[frag_name]:
                            rdf_map[frag_name][end_no] = 0
                        rdf_map[frag_name][end_no] += 1

                elif alignment_allele == variant.alt:
                    counts[CountType.AD] += 1
                    if not aln.is_reverse:
                        counts[CountType.ADP] += 1

                    if self.config.output_fragment_count and frag_name is not None:
                        if frag_name not in adf_map:
                            adf_map[frag_name] = {}
                        if end_no not in adf_map[frag_name]:
                            adf_map[frag_name][end_no] = 0
                        adf_map[frag_name][end_no] += 1

        # Calculate fragment counts
        if self.config.output_fragment_count:
            counts[CountType.DPF] = len(dpf_map)

            fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 0
            fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 0

            for frag_name, end_counts in dpf_map.items():
                # Check for overlapping multimapped reads
                overlap_multimap = False
                for count in end_counts.values():
                    if count > 1:
                        if (
                            self.warning_counts["overlapping_multimap"]
                            < self.config.max_warning_per_type
                        ):
                            logger.warning(
                                f"Fragment {frag_name} has overlapping multiple mapped alignment "
                                f"at site: {variant.chrom}:{variant.pos}"
                            )
                            self.warning_counts["overlapping_multimap"] += 1
                        overlap_multimap = True
                        break

                if overlap_multimap:
                    continue

                # Count fragment ref/alt
                has_ref = frag_name in rdf_map
                has_alt = frag_name in adf_map

                if has_ref and has_alt:
                    # Both ref and alt in fragment
                    counts[CountType.RDF] += fragment_ref_weight
                    counts[CountType.ADF] += fragment_alt_weight
                elif has_ref:
                    counts[CountType.RDF] += 1
                elif has_alt:
                    counts[CountType.ADF] += 1

    def count_bases_snp_numba(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count SNP bases using fast numba_counter for simple SNPs.

        This method converts pysam alignments to NumPy arrays and uses
        the optimized numba_counter for maximum speed on simple SNPs.

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts
        """
        try:
            from .numba_counter import count_snp_base

            # Convert alignments to NumPy arrays for numba_counter
            query_bases = []
            query_qualities = []
            reference_positions = []
            is_reverse_flags = []

            for aln in alignments:
                # Find the base at variant position
                for read_idx, ref_idx in aln.get_aligned_pairs(matches_only=False):
                    if ref_idx == variant.pos:
                        if (
                            aln.query_sequence is not None
                            and aln.query_qualities is not None
                            and read_idx is not None
                        ):
                            query_bases.append(aln.query_sequence[read_idx])
                            query_qualities.append(aln.query_qualities[read_idx])
                            reference_positions.append(ref_idx)
                            is_reverse_flags.append(aln.is_reverse)
                        break

            if not query_bases:
                # No valid bases found, fallback to regular counting
                raise ValueError("No valid bases found at variant position")

            # Convert to NumPy arrays
            import numpy as np

            bases_array = np.array(query_bases, dtype="U1")  # Unicode string array
            quals_array = np.array(query_qualities, dtype=np.uint8)
            pos_array = np.array(reference_positions, dtype=np.int32)
            reverse_array = np.array(is_reverse_flags, dtype=bool)

            # Use numba_counter for fast counting
            dp, rd, ad, dpp, rdp, adp = count_snp_base(
                bases_array,
                quals_array,
                pos_array,
                reverse_array,
                variant.pos,
                variant.ref,
                variant.alt,
                self.config.base_quality_threshold,
            )

            # Convert to our count format
            counts = np.zeros(len(CountType), dtype=np.float32)
            counts[CountType.DP] = dp
            counts[CountType.RD] = rd
            counts[CountType.AD] = ad
            counts[CountType.DPP] = dpp
            counts[CountType.RDP] = rdp
            counts[CountType.ADP] = adp

            # Handle fragment counting if enabled
            if self.config.output_fragment_count:
                # Calculate fragment counts from alignments
                fragment_counts = self._calculate_fragment_counts_numba(alignments, variant)
                counts[CountType.DPF] = fragment_counts["dpf"]
                counts[CountType.RDF] = fragment_counts["rdf"]
                counts[CountType.ADF] = fragment_counts["adf"]

            # Calculate strand bias for this sample (on-the-fly during output)
            # Note: We don't store strand bias in variant object anymore, calculate during output

            # Calculate fragment strand bias if fragment counting is enabled
            if self.config.output_fragment_count:
                # Note: Fragment strand bias uses same forward/reverse as normal strand bias
                # since fragments inherit strand orientation from reads
                pass
            if sample_name not in variant.base_count:
                variant.base_count[sample_name] = counts
            else:
                variant.base_count[sample_name] += counts

        except ImportError:
            # numba not available, fallback to regular counting
            logger.warning("numba_counter not available, falling back to regular SNP counting")
            self.count_bases_snp(variant, alignments, sample_name)
        except Exception as e:
            # Any other error, fallback to regular counting
            logger.warning(f"numba_counter failed: {e}, falling back to regular SNP counting")
            self.count_bases_snp(variant, alignments, sample_name)

    def _calculate_fragment_counts_numba(
        self, alignments: list[pysam.AlignedSegment], variant: VariantEntry
    ) -> dict:
        """
        Calculate fragment counts for numba_counter results.

        Args:
            alignments: Alignments for fragment counting
            variant: Variant for context

        Returns:
            Dictionary with fragment counts
        """
        dpf_map: dict[str, dict[int, int]] = {}
        rdf_map: dict[str, dict[int, int]] = {}
        adf_map: dict[str, dict[int, int]] = {}

        # Group by fragment and track ref/alt
        for aln in alignments:
            if aln.query_name is None:
                continue

            # Find base at variant position for this alignment
            base = None
            for read_idx, ref_idx in aln.get_aligned_pairs(matches_only=False):
                if ref_idx == variant.pos and read_idx is not None:
                    if aln.query_sequence is not None:
                        base = aln.query_sequence[read_idx].upper()
                    break

            if base is None:
                continue

            end_no = 1 if aln.is_read1 else 2
            frag_name = aln.query_name

            # Track depth per fragment
            if frag_name not in dpf_map:
                dpf_map[frag_name] = {}
            dpf_map[frag_name][end_no] = dpf_map[frag_name].get(end_no, 0) + 1

            # Track ref/alt per fragment
            if base == variant.ref:
                if frag_name not in rdf_map:
                    rdf_map[frag_name] = {}
                rdf_map[frag_name][end_no] = rdf_map[frag_name].get(end_no, 0) + 1
            elif base == variant.alt:
                if frag_name not in adf_map:
                    adf_map[frag_name] = {}
                adf_map[frag_name][end_no] = adf_map[frag_name].get(end_no, 0) + 1

        # Calculate final fragment counts
        dpf = len(dpf_map)

        fragment_ref_weight = 0.5 if self.config.fragment_fractional_weight else 1.0
        fragment_alt_weight = 0.5 if self.config.fragment_fractional_weight else 1.0

        rdf = 0.0
        adf = 0.0

        for frag_name, end_counts in dpf_map.items():
            # Check for overlapping multimapped reads
            if any(count > 1 for count in end_counts.values()):
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

        return {"dpf": dpf, "rdf": rdf, "adf": adf}

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

        Args:
            ref_forward: Reference allele count on forward strand
            ref_reverse: Reference allele count on reverse strand
            alt_forward: Alternate allele count on forward strand
            alt_reverse: Alternate allele count on reverse strand
            min_depth: Minimum total depth to calculate bias

        Returns:
            Tuple of (p_value, odds_ratio, bias_direction)
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

    def get_strand_counts_for_sample(
        self, variant: VariantEntry, sample_name: str
    ) -> tuple[int, int, int, int]:
        """
        Get strand-specific counts for a sample.

        Args:
            variant: Variant with counts
            sample_name: Sample name

        Returns:
            Tuple of (ref_forward, ref_reverse, alt_forward, alt_reverse)
        """
        # Get strand-specific counts
        ref_forward = int(variant.get_count(sample_name, CountType.RDP))
        alt_forward = int(variant.get_count(sample_name, CountType.ADP))

        # Calculate reverse strand counts (total - forward)
        ref_total = int(variant.get_count(sample_name, CountType.RD))
        alt_total = int(variant.get_count(sample_name, CountType.AD))

        ref_reverse = max(0, ref_total - ref_forward)
        alt_reverse = max(0, alt_total - alt_forward)

        return ref_forward, ref_reverse, alt_forward, alt_reverse

    def is_simple_snp(self, variant: VariantEntry) -> bool:
        """
        Determine if variant is a simple SNP that can use fast numba_counter.

        Simple SNPs are single-base substitutions that don't require complex
        CIGAR parsing for accurate counting.

        Args:
            variant: Variant to analyze

        Returns:
            True if variant can use numba_counter, False if needs counter.py
        """
        # Must be a SNP (not insertion, deletion, or DNP)
        if not variant.snp:
            return False

        # Must have single-base ref and alt alleles
        if len(variant.ref) != 1 or len(variant.alt) != 1:
            return False

        # Check if we have alignments to analyze
        if not hasattr(variant, "alignments") or not variant.alignments:
            return False

        # For simple SNPs, numba_counter should work well
        # We'll validate this assumption with position checking
        return True

    def validate_variant_position(self, variant: VariantEntry, alignments: list) -> bool:
        """
        Validate that variant position makes sense for counting.

        Args:
            variant: Variant to validate
            alignments: Alignments overlapping variant

        Returns:
            True if position is valid for counting
        """
        if not alignments:
            return False

        # Check if any alignment actually overlaps the variant position
        variant_covered = False
        for aln in alignments:
            if aln.reference_end is not None and aln.reference_start is not None:
                if aln.reference_start <= variant.pos <= aln.reference_end:
                    variant_covered = True
                    break

        return variant_covered

    def smart_count_variant(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Smart counting strategy that chooses optimal algorithm based on variant complexity.

        Strategy:
        - Simple SNPs: Use numba_counter for speed (50-100x faster)
        - Complex variants: Use counter.py for accuracy (better CIGAR handling)

        Args:
            variant: Variant entry to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name for storing counts
        """
        # Validate variant position first
        if not self.validate_variant_position(variant, alignments):
            logger.warning(
                f"No alignments cover variant position for {variant.chrom}:{variant.pos + 1}"
            )
            return

        # Choose counting strategy based on variant complexity and user preference
        if self.config.generic_counting:
            # User explicitly requested generic counting for all variants
            self.count_bases_generic(variant, alignments, sample_name)
        elif self.is_simple_snp(variant):
            try:
                # Try numba_counter for simple SNPs
                self.count_bases_snp_numba(variant, alignments, sample_name)
            except Exception as e:
                logger.warning(
                    f"numba_counter failed for {variant.chrom}:{variant.pos + 1}, falling back to counter.py: {e}"
                )
                # Fallback to counter.py
                self.count_bases_snp(variant, alignments, sample_name)
        else:
            # Use counter.py for complex variants
            if variant.dnp:
                self.count_bases_dnp(variant, alignments, sample_name)
            elif variant.insertion or variant.deletion:
                self.count_bases_indel(variant, alignments, sample_name)
            else:
                # Default to generic counting for unknown variant types
                self.count_bases_generic(variant, alignments, sample_name)

    def count_variant(
        self, variant: VariantEntry, alignments: list[pysam.AlignedSegment], sample_name: str
    ) -> None:
        """
        Count bases for a variant (dispatches to appropriate method).

        Now uses smart counting strategy by default.

        Args:
            variant: Variant to count
            alignments: List of alignments overlapping the variant
            sample_name: Sample name
        """
        # Use smart counting strategy for optimal performance/accuracy balance
        self.smart_count_variant(variant, alignments, sample_name)
