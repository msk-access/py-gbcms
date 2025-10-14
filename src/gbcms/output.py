import logging
from abc import ABC, abstractmethod
from datetime import datetime

from .config import CountType
from .counter import BaseCounter
from .template_engine import template_engine
from .variant import VariantEntry

logger = logging.getLogger(__name__)


"""Output formatting for variant counts."""


class OutputWriter(ABC):
    """Abstract base class for all output format writers."""

    def __init__(self, config, sample_order: list[str]):
        """
        Initialize output writer.

        Args:
            config: Configuration object
            sample_order: Ordered list of sample names
        """
        self.config = config
        self.sample_order = sample_order
        self.template_engine = template_engine

    @abstractmethod
    def write_header(self, output_file: str) -> None:
        """Write format-specific header to output file."""
        pass

    @abstractmethod
    @abstractmethod
    def write_variant(
        self, variant: VariantEntry, output_file: str | None, write_headers: bool = True
    ) -> None:
        pass

    def write_variants(self, variants: list[VariantEntry], output_file: str | None = None) -> None:
        self.validate_variants(variants)
        logger.info(
            f"Writing {len(variants)} variants to {self.get_format_name()} output: {output_file}"
        )

        # Default implementation: call write_variant for each
        # Concrete implementations can override for bulk optimization
        for i, variant in enumerate(variants):
            self.write_variant(variant, output_file, write_headers=(i == 0))

    def get_columns(self) -> list[str]:
        """Get columns for this writer's format using templates."""
        format_name = self.get_format_name()
        return self.template_engine.generate_columns(format_name, self.config, self.sample_order)

    def get_format_name(self) -> str:
        """Get the name of this output format."""
        return self.__class__.__name__.replace("Writer", "").lower()

    def validate_variants(self, variants: list[VariantEntry]) -> None:
        """Validate that variants are compatible with this output format."""
        if not variants:
            raise ValueError(f"No variants provided for {self.get_format_name()} output")

        for variant in variants:
            if not isinstance(variant, VariantEntry):
                raise ValueError(f"Invalid variant type for {self.get_format_name()} output")


class VCFWriter(OutputWriter):
    """VCF format output writer."""

    def __init__(self, config, sample_order: list[str]):
        """Initialize VCF writer."""
        super().__init__(config, sample_order)
        self.counter = BaseCounter(config)

    def write_header(self, output_file: str) -> None:
        """Write VCF header to output file."""
        logger.info(f"Writing VCF header to: {output_file}")

        write_headers = True
        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write("##source=py-gbcms\n")
            f.write(
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">\n'
            )
            f.write(
                '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DP_FORWARD,Number=1,Type=Integer,Description="Forward strand depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RD_FORWARD,Number=1,Type=Integer,Description="Forward strand reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD_FORWARD,Number=1,Type=Integer,Description="Forward strand alternate depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DP_REVERSE,Number=1,Type=Integer,Description="Reverse strand depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RD_REVERSE,Number=1,Type=Integer,Description="Reverse strand reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD_REVERSE,Number=1,Type=Integer,Description="Reverse strand alternate depth for this sample">\n'
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
                '##FORMAT=<ID=RDF_FORWARD,Number=1,Type=Integer,Description="Forward orientation reference fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RDF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation reference fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADF_FORWARD,Number=1,Type=Integer,Description="Forward orientation alternate fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation alternate fragments for this sample">\n'
            )
            f.write(
                '##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction (aggregated across samples)">\n'
            )
            f.write(
                '##INFO=<ID=FSB,Number=3,Type=Float,Description="Fragment strand bias p-value, odds ratio, direction (when fragment counting enabled)">\n'
            )

            # Write column headers
            header_cols = self.get_columns()
            header_cols.extend(self.sample_order)
            f.write("\t".join(header_cols) + "\n")

    def write_variant(
        self, variant: VariantEntry, output_file: str | None, write_headers: bool = True
    ) -> None:
        """Write a single variant to VCF output file efficiently."""
        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            if write_headers:
                # Write VCF header
                f.write("##fileformat=VCFv4.2\n")
                f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
                f.write("##source=py-gbcms\n")
                f.write(
                    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">\n'
                )
                f.write(
                    '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference allele depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=DP_FORWARD,Number=1,Type=Integer,Description="Forward strand depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=RD_FORWARD,Number=1,Type=Integer,Description="Forward strand reference depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=AD_FORWARD,Number=1,Type=Integer,Description="Forward strand alternate depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=DP_REVERSE,Number=1,Type=Integer,Description="Reverse strand depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=RD_REVERSE,Number=1,Type=Integer,Description="Reverse strand reference depth for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=AD_REVERSE,Number=1,Type=Integer,Description="Reverse strand alternate depth for this sample">\n'
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
                    '##FORMAT=<ID=RDF_FORWARD,Number=1,Type=Integer,Description="Forward orientation reference fragments for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=RDF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation reference fragments for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=ADF_FORWARD,Number=1,Type=Integer,Description="Forward orientation alternate fragments for this sample">\n'
                )
                f.write(
                    '##FORMAT=<ID=ADF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation alternate fragments for this sample">\n'
                )
                f.write(
                    '##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction (aggregated across samples)">\n'
                )
                f.write(
                    '##INFO=<ID=FSB,Number=3,Type=Float,Description="Fragment strand bias p-value, odds ratio, direction (when fragment counting enabled)">\n'
                )

            # Write column headers
            header_cols = self.get_columns()
            header_cols.extend(self.sample_order)
            f.write("\t".join(header_cols) + "\n")

            # Write the single variant (extracted from write_variants logic)
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
                ref_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                ref_reverse = rd - ref_forward
                alt_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                alt_reverse = ad - alt_forward

                sb_pval, sb_or, sb_dir = self.counter.calculate_strand_bias(
                    ref_forward, ref_reverse, alt_forward, alt_reverse
                )
                sample_sb_values.append(f"{sb_pval:.6f}:{sb_or:.3f}:{sb_dir}")

                # Fragment strand bias (if enabled) - use real orientation counts
                if self.config.output_fragment_count:
                    # Use real orientation-specific fragment counts (no artificial distribution)
                    ref_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                    ref_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                    alt_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                    alt_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                    fsb_pval, fsb_or, fsb_dir = self.counter.calculate_strand_bias(
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
                total_sb_or = min(total_sb_or, sb_or) if sb_pval < total_sb_pval else total_sb_or
                total_sb_dir = sb_dir if sb_pval < total_sb_pval else total_sb_dir

            # Build INFO field
            info_parts = [
                f"DP={sum(int(variant.get_count(s, CountType.DP)) for s in self.sample_order)}"
            ]

            if total_sb_pval < 1.0:  # Only include if we have valid strand bias
                info_parts.append(f"SB={total_sb_pval:.6f},{total_sb_or:.3f},{total_sb_dir}")

            if self.config.output_fragment_count and total_fsb_pval < 1.0:
                info_parts.append(f"FSB={total_fsb_pval:.6f},{total_fsb_or:.3f},{total_fsb_dir}")

            info_field = ";".join(info_parts)

            # Build FORMAT field
            format_parts = ["DP", "RD", "AD", "VAF"]

            if self.config.output_strand_count:
                format_parts.extend(["DP_FORWARD", "RD_FORWARD", "AD_FORWARD"])

            if self.config.output_strand_count:
                format_parts.extend(["DP_REVERSE", "RD_REVERSE", "AD_REVERSE"])

            if self.config.output_fragment_count:
                format_parts.extend(
                    [
                        "DPF",
                        "RDF",
                        "ADF",
                        "RDF_FORWARD",
                        "RDF_REVERSE",
                        "ADF_FORWARD",
                        "ADF_REVERSE",
                    ]
                )

            format_parts.extend(["SB"])
            if self.config.output_fragment_count:
                format_parts.extend(["FSB"])

            format_field = ":".join(format_parts)

            # Write variant line
            row = [
                variant.chrom,
                str(variant.original_pos),  # Use original 1-indexed coordinates
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

                # Calculate VAF (Variant Allele Frequency)
                vaf = ad / dp if dp > 0 else 0.0

                sample_data = [str(dp), str(rd), str(ad), f"{vaf:.6f}"]

                if self.config.output_strand_count:
                    dp_forward = int(variant.get_count(sample, CountType.DP_FORWARD))
                    rd_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                    ad_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                    dp_reverse = int(variant.get_count(sample, CountType.DP_REVERSE))
                    rd_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                    ad_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))
                    sample_data.extend(
                        [
                            str(dp_forward),
                            str(rd_forward),
                            str(ad_forward),
                            str(dp_reverse),
                            str(rd_reverse),
                            str(ad_reverse),
                        ]
                    )

                if self.config.output_fragment_count:
                    dpf = int(variant.get_count(sample, CountType.DPF))
                    rdf = int(variant.get_count(sample, CountType.RDF))
                    adf = int(variant.get_count(sample, CountType.ADF))
                    rdf_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                    rdf_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                    adf_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                    adf_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))
                    sample_data.extend(
                        [
                            str(dpf),
                            str(rdf),
                            str(adf),
                            str(rdf_forward),
                            str(rdf_reverse),
                            str(adf_forward),
                            str(adf_reverse),
                        ]
                    )

                # Add strand bias data for this sample
                sample_sb_idx = self.sample_order.index(sample)
                sample_data.append(sample_sb_values[sample_sb_idx])

                if self.config.output_fragment_count:
                    sample_data.append(sample_fsb_values[sample_sb_idx])

                row.append(":".join(sample_data))

    def write_variants(
        self, variants: list[VariantEntry], output_file: str | None | None = None
    ) -> None:
        """Write multiple variants to VCF output file."""
        self.validate_variants(variants)
        logger.info(f"Writing {len(variants)} variants to VCF output: {output_file}")

        # Write headers for the first variant only
        write_headers = True
        mode = "w" if write_headers else "a"
        with open(output_file or self.config.output_file, mode) as f:
            # Write header directly to file handle
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write("##source=py-gbcms\n")
            f.write(
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">\n'
            )
            f.write(
                '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DP_FORWARD,Number=1,Type=Integer,Description="Forward strand depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RD_FORWARD,Number=1,Type=Integer,Description="Forward strand reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD_FORWARD,Number=1,Type=Integer,Description="Forward strand alternate depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=DP_REVERSE,Number=1,Type=Integer,Description="Reverse strand depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RD_REVERSE,Number=1,Type=Integer,Description="Reverse strand reference depth for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=AD_REVERSE,Number=1,Type=Integer,Description="Reverse strand alternate depth for this sample">\n'
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
                '##FORMAT=<ID=RDF_FORWARD,Number=1,Type=Integer,Description="Forward orientation reference fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=RDF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation reference fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADF_FORWARD,Number=1,Type=Integer,Description="Forward orientation alternate fragments for this sample">\n'
            )
            f.write(
                '##FORMAT=<ID=ADF_REVERSE,Number=1,Type=Integer,Description="Reverse orientation alternate fragments for this sample">\n'
            )
            f.write(
                '##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction (aggregated across samples)">\n'
            )
            f.write(
                '##INFO=<ID=FSB,Number=3,Type=Float,Description="Fragment strand bias p-value, odds ratio, direction (when fragment counting enabled)">\n'
            )

            # Write column headers
            header_cols = self.get_columns()
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
                    ref_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                    ref_reverse = rd - ref_forward
                    alt_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                    alt_reverse = ad - alt_forward

                    sb_pval, sb_or, sb_dir = self.counter.calculate_strand_bias(
                        ref_forward, ref_reverse, alt_forward, alt_reverse
                    )
                    sample_sb_values.append(f"{sb_pval:.6f}:{sb_or:.3f}:{sb_dir}")

                    # Fragment strand bias (if enabled) - use real orientation counts
                    if self.config.output_fragment_count:
                        # Use real orientation-specific fragment counts (no artificial distribution)
                        ref_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        ref_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        alt_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        alt_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                        fsb_pval, fsb_or, fsb_dir = self.counter.calculate_strand_bias(
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
                format_parts = ["DP", "RD", "AD", "VAF"]

                if self.config.output_strand_count:
                    format_parts.extend(["DP_FORWARD", "RD_FORWARD", "AD_FORWARD"])

                if self.config.output_strand_count:
                    format_parts.extend(["DP_REVERSE", "RD_REVERSE", "AD_REVERSE"])

                if self.config.output_fragment_count:
                    format_parts.extend(
                        [
                            "DPF",
                            "RDF",
                            "ADF",
                            "RDF_FORWARD",
                            "RDF_REVERSE",
                            "ADF_FORWARD",
                            "ADF_REVERSE",
                        ]
                    )

                format_parts.extend(["SB"])
                if self.config.output_fragment_count:
                    format_parts.extend(["FSB"])

                format_field = ":".join(format_parts)

                # Write variant line
                row = [
                    variant.chrom,
                    str(variant.original_pos),  # Use original 1-indexed coordinates
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

                    # Calculate VAF (Variant Allele Frequency)
                    vaf = ad / dp if dp > 0 else 0.0

                    sample_data = [str(dp), str(rd), str(ad), f"{vaf:.6f}"]

                    if self.config.output_strand_count:
                        dp_forward = int(variant.get_count(sample, CountType.DP_FORWARD))
                        rd_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                        ad_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                        dp_reverse = int(variant.get_count(sample, CountType.DP_REVERSE))
                        rd_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                        ad_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))
                        sample_data.extend(
                            [
                                str(dp_forward),
                                str(rd_forward),
                                str(ad_forward),
                                str(dp_reverse),
                                str(rd_reverse),
                                str(ad_reverse),
                            ]
                        )

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        rdf_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        rdf_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        adf_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        adf_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))
                        sample_data.extend(
                            [
                                str(dpf),
                                str(rdf),
                                str(adf),
                                str(rdf_forward),
                                str(rdf_reverse),
                                str(adf_forward),
                                str(adf_reverse),
                            ]
                        )

                    # Add strand bias data for this sample
                    sample_sb_idx = self.sample_order.index(sample)
                    sample_data.append(sample_sb_values[sample_sb_idx])

                    if self.config.output_fragment_count:
                        sample_data.append(sample_fsb_values[sample_sb_idx])

                    row.append(":".join(sample_data))

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to VCF output file")


class SampleAgnosticMAFWriter(OutputWriter):
    """Sample-agnostic MAF format output writer (one row per sample per variant)."""

    def __init__(self, config, sample_order: list[str]):
        """Initialize Sample-Agnostic MAF writer."""
        super().__init__(config, sample_order)
        self.counter = BaseCounter(config)

    def write_header(self, output_file: str) -> None:
        """Write MAF header to output file."""
        logger.info(f"Writing sample-agnostic MAF header to: {output_file}")

        write_headers = True
        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            # Write header - standard MAF columns plus tumor-only count columns
            columns = self.get_columns()
            header_cols = columns
            f.write("\t".join(header_cols) + "\n")

    def write_variant(
        self, variant: VariantEntry, output_file: str | None, write_headers: bool = True
    ) -> None:
        """Write a single variant to MAF output file."""
        # For single variant, we need to write header + variant for all samples
        # This is less efficient but maintains API compatibility
        variants = [variant]
        self.write_variants(variants)

    def write_variants(
        self, variants: list[VariantEntry], output_file: str | None | None = None
    ) -> None:
        """Write multiple variants to sample-agnostic MAF output file."""
        self.validate_variants(variants)
        logger.info(
            f"Writing {len(variants)} variants to sample-agnostic MAF output: {self.config.output_file}"
        )

        # Write headers for the first variant only
        write_headers = True
        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            # Write header - standard MAF columns plus tumor-only count columns
            columns = self.get_columns()
            columns = self.get_columns()
            f.write("\t".join(columns) + "\n")

            # Write variants - one row per sample per variant
            for variant in variants:
                # Base MAF columns (same for all samples)
                base_cols = [
                    variant.gene,
                    variant.chrom,
                    str(variant.original_pos),  # Use original 1-indexed coordinates
                    str(variant.original_end_pos),
                    variant.maf_ref,
                    variant.maf_ref,  # Tumor_Seq_Allele1 = reference
                    variant.maf_alt if variant.maf_alt else "",  # Tumor_Seq_Allele2 = alt (variant)
                    "",  # Tumor_Sample_Barcode (will be set per sample)
                    "",  # Matched_Norm_Sample_Barcode (always empty)
                    variant.effect,
                ]

                # Generate one row per sample
                for sample in self.sample_order:
                    row = base_cols.copy()
                    row[7] = sample  # Set Tumor_Sample_Barcode to sample name

                    # Get counts for this sample
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))
                    vaf = ad / dp if dp > 0 else 0.0

                    # Add tumor counts (sample-specific)
                    sample_counts = [
                        str(dp),  # t_depth
                        str(rd),  # t_ref_count
                        str(ad),  # t_alt_count
                        f"{vaf:.6f}",  # t_vaf
                    ]

                    if self.config.output_strand_count:
                        dp_forward = int(variant.get_count(sample, CountType.DP_FORWARD))
                        rd_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                        ad_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                        dp_reverse = int(variant.get_count(sample, CountType.DP_REVERSE))
                        rd_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                        ad_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))

                        sample_counts.extend(
                            [
                                str(dp_forward),  # t_depth_forward
                                str(rd_forward),  # t_ref_count_forward
                                str(ad_forward),  # t_alt_count_forward
                                str(dp_reverse),  # t_depth_reverse
                                str(rd_reverse),  # t_ref_count_reverse
                                str(ad_reverse),  # t_alt_count_reverse
                            ]
                        )

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        rdf_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        rdf_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        adf_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        adf_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                        sample_counts.extend(
                            [
                                str(dpf),  # t_depth_fragment
                                str(rdf),  # t_ref_count_fragment
                                str(adf),  # t_alt_count_fragment
                                str(rdf_forward),  # t_ref_count_fragment_forward
                                str(rdf_reverse),  # t_ref_count_fragment_reverse
                                str(adf_forward),  # t_alt_count_fragment_forward
                                str(adf_reverse),  # t_alt_count_fragment_reverse
                            ]
                        )

                    # Calculate and add strand bias for this sample (separate columns)
                    ref_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                    ref_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                    alt_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                    alt_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))

                    sb_pval, sb_or, sb_dir = self.counter.calculate_strand_bias(
                        ref_forward, ref_reverse, alt_forward, alt_reverse
                    )

                    sample_counts.extend(
                        [
                            f"{sb_pval:.6f}",  # t_strand_bias_pval
                            f"{sb_or:.3f}",  # t_strand_bias_or
                            sb_dir,  # t_strand_bias_dir
                        ]
                    )

                    if self.config.output_fragment_count:
                        # Calculate fragment strand bias
                        ref_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        ref_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        alt_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        alt_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                        fsb_pval, fsb_or, fsb_dir = self.counter.calculate_strand_bias(
                            ref_forward, ref_reverse, alt_forward, alt_reverse
                        )

                        sample_counts.extend(
                            [
                                f"{fsb_pval:.6f}",  # t_fragment_strand_bias_pval
                                f"{fsb_or:.3f}",  # t_fragment_strand_bias_or
                                fsb_dir,  # t_fragment_strand_bias_dir
                            ]
                        )

                    row.extend(sample_counts)
                    f.write("\t".join(row) + "\n")

        logger.info(
            f"Successfully wrote {len(variants) * len(self.sample_order)} variant-sample combinations to MAF output file"
        )


class FilloutWriter(OutputWriter):
    """Fillout MAF format output writer (one row per variant with all samples)."""

    def __init__(self, config, sample_order: list[str]):
        """Initialize Fillout writer."""
        super().__init__(config, sample_order)
        self.counter = BaseCounter(config)

    def write_header(self, output_file: str) -> None:
        """Write fillout MAF header to output file."""
        logger.info(f"Writing fillout MAF header to: {output_file}")

        write_headers = True
        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            # Write header - standard MAF columns plus sample-specific count columns
            columns = self.get_columns()
            header_cols = columns

            # Add count columns for each sample
            for _sample in self.sample_order:
                columns = self.get_columns()
                f.write("\t".join(header_cols + ["\t".join(columns)]) + "\n")

    def write_variant(
        self, variant: VariantEntry, output_file: str | None, write_headers: bool = True
    ) -> None:
        """Write a single variant to fillout MAF output file."""
        # For single variant, we need to write header + variant
        # This is less efficient but maintains API compatibility
        variants = [variant]
        self.write_variants(variants)

    def write_variants(
        self, variants: list[VariantEntry], output_file: str | None | None = None
    ) -> None:
        """Write multiple variants to fillout MAF output file."""
        self.validate_variants(variants)
        logger.info(f"Writing {len(variants)} variants to fillout MAF output: {output_file}")
        write_headers = True

        mode = "w" if write_headers else "a"
        with open(self.config.output_file, mode) as f:
            # Write header - standard MAF columns plus sample-specific count columns
            columns = self.get_columns()
            # Add count columns for each sample

            # Add count columns for each sample
            for _sample in self.sample_order:
                columns = self.get_columns()
                f.write("\t".join(columns) + "\n")

            # Write variants - one row per variant with all samples
            for variant in variants:
                row = [
                    variant.gene,
                    variant.chrom,
                    str(variant.original_pos),  # Use original 1-indexed coordinates
                    str(variant.original_end_pos),
                    variant.maf_ref,
                    variant.maf_ref,  # Tumor_Seq_Allele1 = reference
                    variant.maf_alt if variant.maf_alt else "",  # Tumor_Seq_Allele2 = alt (variant)
                    "",  # Tumor_Sample_Barcode (empty for fillout format)
                    "",  # Matched_Norm_Sample_Barcode (empty for fillout format)
                    variant.effect,
                ]

                # Add counts for each sample in this variant's row
                for sample in self.sample_order:
                    dp = int(variant.get_count(sample, CountType.DP))
                    rd = int(variant.get_count(sample, CountType.RD))
                    ad = int(variant.get_count(sample, CountType.AD))
                    vaf = ad / dp if dp > 0 else 0.0

                    # Add basic counts for this sample
                    row.extend([str(dp), str(rd), str(ad), f"{vaf:.6f}"])

                    if self.config.output_strand_count:
                        dp_forward = int(variant.get_count(sample, CountType.DP_FORWARD))
                        rd_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                        ad_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                        dp_reverse = int(variant.get_count(sample, CountType.DP_REVERSE))
                        rd_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                        ad_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))

                        row.extend(
                            [
                                str(dp_forward),
                                str(rd_forward),
                                str(ad_forward),
                                str(dp_reverse),
                                str(rd_reverse),
                                str(ad_reverse),
                            ]
                        )

                    if self.config.output_fragment_count:
                        dpf = int(variant.get_count(sample, CountType.DPF))
                        rdf = int(variant.get_count(sample, CountType.RDF))
                        adf = int(variant.get_count(sample, CountType.ADF))
                        rdf_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        rdf_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        adf_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        adf_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                        row.extend(
                            [
                                str(dpf),
                                str(rdf),
                                str(adf),
                                str(rdf_forward),
                                str(rdf_reverse),
                                str(adf_forward),
                                str(adf_reverse),
                            ]
                        )

                    # Calculate and add strand bias for this sample (separate columns)
                    ref_forward = int(variant.get_count(sample, CountType.RD_FORWARD))
                    ref_reverse = int(variant.get_count(sample, CountType.RD_REVERSE))
                    alt_forward = int(variant.get_count(sample, CountType.AD_FORWARD))
                    alt_reverse = int(variant.get_count(sample, CountType.AD_REVERSE))

                    sb_pval, sb_or, sb_dir = self.counter.calculate_strand_bias(
                        ref_forward, ref_reverse, alt_forward, alt_reverse
                    )

                    row.extend([f"{sb_pval:.6f}", f"{sb_or:.3f}", sb_dir])

                    if self.config.output_fragment_count:
                        # Calculate fragment strand bias
                        ref_forward = int(variant.get_count(sample, CountType.RDF_FORWARD))
                        ref_reverse = int(variant.get_count(sample, CountType.RDF_REVERSE))
                        alt_forward = int(variant.get_count(sample, CountType.ADF_FORWARD))
                        alt_reverse = int(variant.get_count(sample, CountType.ADF_REVERSE))

                        fsb_pval, fsb_or, fsb_dir = self.counter.calculate_strand_bias(
                            ref_forward, ref_reverse, alt_forward, alt_reverse
                        )

                        row.extend([f"{fsb_pval:.6f}", f"{fsb_or:.3f}", fsb_dir])

                f.write("\t".join(row) + "\n")

        logger.info(f"Successfully wrote {len(variants)} variants to fillout output file")
