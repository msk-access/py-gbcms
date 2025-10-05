"""Output formatting for variant counts."""

import logging
from datetime import datetime

from .config import Config, CountType
from .variant import VariantEntry

logger = logging.getLogger(__name__)


class OutputFormatter:
    """Formats and writes output files."""

    def __init__(self, config: Config, sample_order: list[str]):
        """
        Initialize output formatter.

        Args:
            config: Configuration object
            sample_order: Ordered list of sample names
        """
        self.config = config
        self.sample_order = sample_order

    def write_vcf_output(self, variants: list[VariantEntry]) -> None:
        """
        Write output in proper VCF format with strand bias in INFO field.

        Args:
            variants: List of variants with counts and strand bias
        """
        logger.info(f"Writing VCF output to: {self.config.output_file}")

        with open(self.config.output_file, "w") as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write("##source=py-gbcms\n")
            f.write(
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">\n'
            )
            f.write(
                '##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction (aggregated across samples)">\n'
            )
            f.write(
                '##INFO=<ID=FSB,Number=3,Type=Float,Description="Fragment strand bias p-value, odds ratio, direction (when fragment counting enabled)">\n'
            )
            f.write(
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DPP,Number=1,Type=Integer,Description="Positive strand depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RDP,Number=1,Type=Integer,Description="Positive strand reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADP,Number=1,Type=Integer,Description="Positive strand alternate depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Fragment depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Fragment reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Fragment alternate depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=FSB,Number=3,Type=Float,Description="Fragment strand bias p-value, odds ratio, direction for this sample">\n'
            )

            # Write column headers
            header_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            header_cols.extend(self.sample_order)
            f.write("\t".join(header_cols) + "\n")

            # Write variants
            for variant in variants:
                # Calculate aggregate strand bias across all samples for INFO field
                total_sb_pval = 1.0
                total_sb_or = 1.0
                total_sb_dir = "none"
                total_fsb_pval = 1.0
                total_fsb_or = 1.0
                total_fsb_dir = "none"

                sample_sb_values = []
                sample_fsb_values = []

                for sample in self.sample_order:
                    # Get counts
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))

                    # Calculate strand bias for this sample
                    ref_forward = int(variant.get_count(sample, CountType.RDP))
                    ref_reverse = rd - ref_forward
                    alt_forward = int(variant.get_count(sample, CountType.ADP))
                    alt_reverse = ad - alt_forward

                    sb_pval, sb_or, sb_dir = self._calculate_strand_bias_for_output(
                        ref_forward, ref_reverse, alt_forward, alt_reverse
                    )
                    sample_sb_values.append(f"{sb_pval:.6f}:{sb_or:.3f}:{sb_dir}")

                    # Fragment strand bias (if enabled)
                    if self.config.output_fragment_count:
                        fsb_pval, fsb_or, fsb_dir = self._calculate_strand_bias_for_output(
                            ref_forward, ref_reverse, alt_forward, alt_reverse
                        )
                        sample_fsb_values.append(f"{fsb_pval:.6f}:{fsb_or:.3f}:{fsb_dir}")
                        total_fsb_pval = min(total_fsb_pval, fsb_pval)
                        total_fsb_or = (
                            min(total_fsb_or, fsb_or) if fsb_pval < total_fsb_pval else total_fsb_or
                        )
                        total_fsb_dir = fsb_dir if fsb_pval < total_fsb_pval else total_fsb_dir
                    else:
                        sample_fsb_values.append(".:.:none")

                    # Update aggregate values (use minimum p-value as most significant)
                    total_sb_pval = min(total_sb_pval, sb_pval)
                    total_sb_or = (
                        min(total_sb_or, sb_or) if sb_pval < total_sb_pval else total_sb_or
                    )
                    total_sb_dir = sb_dir if sb_pval < total_sb_pval else total_sb_dir

                # Build INFO field
                info_parts = [
                    f"DP={sum(int(variant.get_count(s, CountType.DP)) for s in self.sample_order)}"
                ]

                if total_sb_pval < 1.0:  # Only include if we have valid strand bias
                    info_parts.append(f"SB={total_sb_pval:.6f},{total_sb_or:.3f},{total_sb_dir}")

                if self.config.output_fragment_count and total_fsb_pval < 1.0:
                    info_parts.append(
                        f"FSB={total_fsb_pval:.6f},{total_fsb_or:.3f},{total_fsb_dir}"
                    )

                info_field = ";".join(info_parts)

                # Build FORMAT field
                format_parts = ["DP", "RD", "AD"]

                if self.config.output_positive_count:
                    format_parts.extend(["DPP", "RDP", "ADP"])

                if self.config.output_fragment_count:
                    format_parts.extend(["DPF", "RDF", "ADF"])

                format_parts.extend(["SB"])
                if self.config.output_fragment_count:
                    format_parts.extend(["FSB"])

                format_field = ":".join(format_parts)

                # Write variant line
                row = [
                    variant.chrom,
                    str(variant.pos + 1),  # Convert to 1-indexed
                    ".",  # ID
                    variant.ref,
                    variant.alt,
                    ".",  # QUAL
                    ".",  # FILTER
                    info_field,
                    format_field,
                ]

                # Add sample data
                for sample in self.sample_order:
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))

                    sample_data = [str(dp), str(rd), str(ad)]

                    if self.config.output_positive_count:
                        dpp = int(variant.get_count(sample, CountType.DPP))
                        rdp = int(variant.get_count(sample, CountType.RDP))
                        adp = int(variant.get_count(sample, CountType.ADP))
                        sample_data.extend([str(dpp), str(rdp), str(adp)])

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        sample_data.extend([str(dpf), str(rdf), str(adf)])

                    # Add strand bias data for this sample
                    sample_sb_idx = self.sample_order.index(sample)
                    sample_data.append(sample_sb_values[sample_sb_idx])

                    if self.config.output_fragment_count:
                        sample_data.append(sample_fsb_values[sample_sb_idx])

                    row.append(":".join(sample_data))

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to VCF output file")

    def _calculate_strand_bias_for_output(
        self,
        ref_forward: int,
        ref_reverse: int,
        alt_forward: int,
        alt_reverse: int,
        min_depth: int = 10,
    ) -> tuple[float, float, str]:
        """
        Calculate strand bias using Fisher's exact test for output.

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

    def write_maf_output(self, variants: list[VariantEntry]) -> None:
        """
        Write output in MAF format.

        Args:
            variants: List of variants with counts
        """
        logger.info(f"Writing MAF output to: {self.config.output_file}")

        with open(self.config.output_file, "w") as f:
            # Write header (use first variant's MAF line to get column names)
            if variants and variants[0].maf_line:
                # Parse header from first MAF line structure
                header_cols = [
                    "Hugo_Symbol",
                    "Chromosome",
                    "Start_Position",
                    "End_Position",
                    "Reference_Allele",
                    "Tumor_Seq_Allele1",
                    "Tumor_Seq_Allele2",
                    "Tumor_Sample_Barcode",
                    "Matched_Norm_Sample_Barcode",
                    "Variant_Classification",
                ]

                # Add count columns for tumor and normal
                count_cols = [
                    "t_depth",
                    "t_ref_count",
                    "t_alt_count",
                    "n_depth",
                    "n_ref_count",
                    "n_alt_count",
                ]

                if self.config.output_positive_count:
                    count_cols.extend(
                        [
                            "t_depth_forward",
                            "t_ref_count_forward",
                            "t_alt_count_forward",
                            "n_depth_forward",
                            "n_ref_count_forward",
                            "n_alt_count_forward",
                        ]
                    )

                if self.config.output_fragment_count:
                    count_cols.extend(
                        [
                            "t_depth_fragment",
                            "t_ref_count_fragment",
                            "t_alt_count_fragment",
                            "n_depth_fragment",
                            "n_ref_count_fragment",
                            "n_alt_count_fragment",
                        ]
                    )

                # Add strand bias columns for tumor and normal samples
                count_cols.extend(
                    [
                        "t_strand_bias_pval",
                        "t_strand_bias_or",
                        "t_strand_bias_dir",
                        "n_strand_bias_pval",
                        "n_strand_bias_or",
                        "n_strand_bias_dir",
                    ]
                )

                if self.config.output_fragment_count:
                    count_cols.extend(
                        [
                            "t_fragment_strand_bias_pval",
                            "t_fragment_strand_bias_or",
                            "t_fragment_strand_bias_dir",
                            "n_fragment_strand_bias_pval",
                            "n_fragment_strand_bias_or",
                            "n_fragment_strand_bias_dir",
                        ]
                    )

                f.write("\t".join(header_cols + count_cols) + "\n")

            # Write variants
            for variant in variants:
                row = [
                    variant.gene,
                    variant.chrom,
                    str(variant.maf_pos + 1),  # Convert back to 1-indexed
                    str(variant.maf_end_pos + 1),
                    variant.maf_ref,
                    variant.maf_alt if variant.maf_alt else "",
                    "",  # Tumor_Seq_Allele2
                    variant.tumor_sample,
                    variant.normal_sample,
                    variant.effect,
                ]

                # Get tumor counts
                t_dp = int(variant.get_count(variant.tumor_sample, CountType.DP))
                t_rd = int(variant.get_count(variant.tumor_sample, CountType.RD))
                t_ad = int(variant.get_count(variant.tumor_sample, CountType.AD))

                # Get normal counts
                n_dp = int(variant.get_count(variant.normal_sample, CountType.DP))
                n_rd = int(variant.get_count(variant.normal_sample, CountType.RD))
                n_ad = int(variant.get_count(variant.normal_sample, CountType.AD))

                row.extend([str(t_dp), str(t_rd), str(t_ad), str(n_dp), str(n_rd), str(n_ad)])

                if self.config.output_positive_count:
                    t_dpp = int(variant.get_count(variant.tumor_sample, CountType.DPP))
                    t_rdp = int(variant.get_count(variant.tumor_sample, CountType.RDP))
                    t_adp = int(variant.get_count(variant.tumor_sample, CountType.ADP))
                    n_dpp = int(variant.get_count(variant.normal_sample, CountType.DPP))
                    n_rdp = int(variant.get_count(variant.normal_sample, CountType.RDP))
                    n_adp = int(variant.get_count(variant.normal_sample, CountType.ADP))
                    row.extend(
                        [str(t_dpp), str(t_rdp), str(t_adp), str(n_dpp), str(n_rdp), str(n_adp)]
                    )

                if self.config.output_fragment_count:
                    t_dpf = int(variant.get_count(variant.tumor_sample, CountType.DPF))
                    t_rdf = int(variant.get_count(variant.tumor_sample, CountType.RDF))
                    t_adf = int(variant.get_count(variant.tumor_sample, CountType.ADF))
                    n_dpf = int(variant.get_count(variant.normal_sample, CountType.DPF))
                    n_rdf = int(variant.get_count(variant.normal_sample, CountType.RDF))
                    n_adf = int(variant.get_count(variant.normal_sample, CountType.ADF))
                    row.extend(
                        [str(t_dpf), str(t_rdf), str(t_adf), str(n_dpf), str(n_rdf), str(n_adf)]
                    )

                # Add strand bias information for tumor and normal
                # Calculate tumor strand bias on-the-fly
                t_ref_forward = int(variant.get_count(variant.tumor_sample, CountType.RDP))
                t_ref_reverse = t_rd - t_ref_forward
                t_alt_forward = int(variant.get_count(variant.tumor_sample, CountType.ADP))
                t_alt_reverse = t_ad - t_alt_forward

                t_sb_pval, t_sb_or, t_sb_dir = self._calculate_strand_bias_for_output(
                    t_ref_forward, t_ref_reverse, t_alt_forward, t_alt_reverse
                )

                # Calculate normal strand bias on-the-fly
                n_ref_forward = int(variant.get_count(variant.normal_sample, CountType.RDP))
                n_ref_reverse = n_rd - n_ref_forward
                n_alt_forward = int(variant.get_count(variant.normal_sample, CountType.ADP))
                n_alt_reverse = n_ad - n_alt_forward

                n_sb_pval, n_sb_or, n_sb_dir = self._calculate_strand_bias_for_output(
                    n_ref_forward, n_ref_reverse, n_alt_forward, n_alt_reverse
                )

                row.extend(
                    [
                        f"{t_sb_pval:.6f}",
                        f"{t_sb_or:.3f}",
                        t_sb_dir,
                        f"{n_sb_pval:.6f}",
                        f"{n_sb_or:.3f}",
                        n_sb_dir,
                    ]
                )

                if self.config.output_fragment_count:
                    # Calculate fragment strand bias for tumor and normal
                    t_fsb_pval, t_fsb_or, t_fsb_dir = self._calculate_strand_bias_for_output(
                        t_ref_forward, t_ref_reverse, t_alt_forward, t_alt_reverse
                    )
                    n_fsb_pval, n_fsb_or, n_fsb_dir = self._calculate_strand_bias_for_output(
                        n_ref_forward, n_ref_reverse, n_alt_forward, n_alt_reverse
                    )

                    row.extend(
                        [
                            f"{t_fsb_pval:.6f}",
                            f"{t_fsb_or:.3f}",
                            t_fsb_dir,
                            f"{n_fsb_pval:.6f}",
                            f"{n_fsb_or:.3f}",
                            n_fsb_dir,
                        ]
                    )

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to MAF output file")

    def write_fillout_output(self, variants: list[VariantEntry]) -> None:
        """
        Write output in fillout format (extended MAF with all samples).

        Args:
            variants: List of variants with counts
        """
        logger.info(f"Writing fillout output to: {self.config.output_file}")

        with open(self.config.output_file, "w") as f:
            # Write header
            header_cols = [
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele1",
                "Tumor_Seq_Allele2",
                "Tumor_Sample_Barcode",
                "Matched_Norm_Sample_Barcode",
                "Variant_Classification",
            ]

            # Add count columns for each sample
            for sample in self.sample_order:
                header_cols.extend([f"{sample}:DP", f"{sample}:RD", f"{sample}:AD"])
                if self.config.output_positive_count:
                    header_cols.extend([f"{sample}:DPP", f"{sample}:RDP", f"{sample}:ADP"])
                if self.config.output_fragment_count:
                    header_cols.extend([f"{sample}:DPF", f"{sample}:RDF", f"{sample}:ADF"])

                # Add strand bias columns for each sample
                header_cols.extend([f"{sample}:SB_PVAL", f"{sample}:SB_OR", f"{sample}:SB_DIR"])

                if self.config.output_fragment_count:
                    header_cols.extend(
                        [f"{sample}:FSB_PVAL", f"{sample}:FSB_OR", f"{sample}:FSB_DIR"]
                    )

            f.write("\t".join(header_cols) + "\n")

            # Write variants
            for variant in variants:
                row = [
                    variant.gene,
                    variant.chrom,
                    str(variant.maf_pos + 1),
                    str(variant.maf_end_pos + 1),
                    variant.maf_ref,
                    variant.maf_alt if variant.maf_alt else "",
                    "",
                    variant.tumor_sample,
                    variant.normal_sample,
                    variant.effect,
                ]

                for sample in self.sample_order:
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))
                    row.extend([str(dp), str(rd), str(ad)])

                    if self.config.output_positive_count:
                        dpp = int(variant.get_count(sample, CountType.DPP))
                        rdp = int(variant.get_count(sample, CountType.RDP))
                        adp = int(variant.get_count(sample, CountType.ADP))
                        row.extend([str(dpp), str(rdp), str(adp)])

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        row.extend([str(dpf), str(rdf), str(adf)])

                    # Add strand bias information for this sample
                    # Calculate strand bias on-the-fly using normal counts
                    sample_ref_forward = int(variant.get_count(sample, CountType.RDP))
                    sample_ref_reverse = rd - sample_ref_forward
                    sample_alt_forward = int(variant.get_count(sample, CountType.ADP))
                    sample_alt_reverse = ad - sample_alt_forward

                    sb_pval, sb_or, sb_dir = self._calculate_strand_bias_for_output(
                        sample_ref_forward,
                        sample_ref_reverse,
                        sample_alt_forward,
                        sample_alt_reverse,
                    )
                    row.extend([f"{sb_pval:.6f}", f"{sb_or:.3f}", sb_dir])

                    if self.config.output_fragment_count:
                        # For fragment strand bias, use the same forward/reverse counts
                        # (fragments inherit strand orientation from their constituent reads)
                        fsb_pval, fsb_or, fsb_dir = self._calculate_strand_bias_for_output(
                            sample_ref_forward,
                            sample_ref_reverse,
                            sample_alt_forward,
                            sample_alt_reverse,
                        )
                        row.extend([f"{fsb_pval:.6f}", f"{fsb_or:.3f}", fsb_dir])

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to fillout output file")
