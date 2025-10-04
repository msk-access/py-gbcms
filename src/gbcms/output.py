"""Output formatting for variant counts."""

import logging

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
        Write output in VCF-like format.

        Args:
            variants: List of variants with counts
        """
        logger.info(f"Writing output to: {self.config.output_file}")

        with open(self.config.output_file, "w") as f:
            # Write header
            header_cols = ["Chrom", "Pos", "Ref", "Alt"]

            for sample in self.sample_order:
                header_cols.extend([f"{sample}:DP", f"{sample}:RD", f"{sample}:AD"])
                if self.config.output_positive_count:
                    header_cols.extend([f"{sample}:DPP", f"{sample}:RDP", f"{sample}:ADP"])
                if self.config.output_negative_count:
                    header_cols.extend([f"{sample}:DPN", f"{sample}:RDN", f"{sample}:ADN"])
                if self.config.output_fragment_count:
                    header_cols.extend([f"{sample}:DPF", f"{sample}:RDF", f"{sample}:ADF"])

            f.write("\t".join(header_cols) + "\n")

            # Write variants
            for variant in variants:
                row = [
                    variant.chrom,
                    str(variant.pos + 1),  # Convert back to 1-indexed
                    variant.ref,
                    variant.alt,
                ]

                for sample in self.sample_order:
                    # Get counts (default to 0 if not present)
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))

                    row.extend([str(dp), str(rd), str(ad)])

                    if self.config.output_positive_count:
                        dpp = int(variant.get_count(sample, CountType.DPP))
                        rdp = int(variant.get_count(sample, CountType.RDP))
                        adp = int(variant.get_count(sample, CountType.ADP))
                        row.extend([str(dpp), str(rdp), str(adp)])

                    if self.config.output_negative_count:
                        # Negative counts = total - positive
                        dpn = dp - int(variant.get_count(sample, CountType.DPP))
                        rdn = rd - int(variant.get_count(sample, CountType.RDP))
                        adn = ad - int(variant.get_count(sample, CountType.ADP))
                        row.extend([str(dpn), str(rdn), str(adn)])

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        row.extend([str(dpf), str(rdf), str(adf)])

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to output file")

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

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to fillout output file")
