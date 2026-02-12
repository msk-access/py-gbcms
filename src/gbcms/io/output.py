"""
Output Writers: Formatting results for VCF and MAF.

This module provides classes to write processed variants and their counts
to output files, handling format-specific columns and headers.
"""

import csv
import logging
from pathlib import Path
from typing import Any

from ..core.kernel import CoordinateKernel
from ..models.core import Variant

__all__ = ["OutputWriter", "MafWriter", "VcfWriter"]

logger = logging.getLogger(__name__)


class OutputWriter:
    """Abstract base class for output writers."""

    def write(self, variant: Variant, counts: Any):
        raise NotImplementedError

    def close(self):
        pass


class MafWriter(OutputWriter):
    """
    Writes results to a MAF-like file (Fillout format).

    Supports two output strategies based on input format:
    - MAF→MAF: Preserves all original MAF columns, appends gbcms count columns.
      Original positional/allele/type columns are NEVER overwritten.
    - VCF→MAF: Generates GDC-compliant MAF coordinates from internal VCF-style
      representation using CoordinateKernel.internal_to_maf().

    Count column names are controlled by the column_prefix parameter:
    - Default (empty): 'ref_count', 'alt_count', 'total_count', etc.
    - Legacy ('t_'):   't_ref_count', 't_alt_count', 't_total_count', etc.
    """

    # Default MAF columns for VCF→MAF output (minimal GDC-compatible set)
    _DEFAULT_MAF_HEADERS = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Strand",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
    ]

    def __init__(self, path: Path, column_prefix: str = "", preserve_barcode: bool = False):
        """
        Initialize MafWriter.

        Args:
            path: Output file path.
            column_prefix: Prefix for gbcms count columns (e.g., '', 't_', 'gbcms_').
            preserve_barcode: If True, keep original Tumor_Sample_Barcode from
                input MAF. If False (default), override with BAM sample name.
                Only applies to MAF→MAF; VCF→MAF always uses BAM name.
        """
        self.path = path
        self.column_prefix = column_prefix
        self.preserve_barcode = preserve_barcode
        self.file = open(path, "w")
        self.writer: csv.DictWriter | None = None
        self._headers_written = False
        logger.debug(
            "MafWriter initialized: path=%s, column_prefix='%s', preserve_barcode=%s",
            path,
            column_prefix,
            preserve_barcode,
        )

    def _gbcms_column_names(self) -> list[str]:
        """
        Build the list of gbcms-generated count column names with the configured prefix.

        Returns:
            Ordered list of gbcms column names.
        """
        p = self.column_prefix
        return [
            # Core counts
            f"{p}ref_count",
            f"{p}alt_count",
            f"{p}total_count",
            f"{p}vaf",
            # Fragment counts
            f"{p}ref_count_fragment",
            f"{p}alt_count_fragment",
            f"{p}total_count_fragment",
            f"{p}vaf_fragment",
            # Strand bias (unprefixed — always unique)
            "strand_bias_p_value",
            "strand_bias_odds_ratio",
            "fragment_strand_bias_p_value",
            "fragment_strand_bias_odds_ratio",
            # Strand counts
            f"{p}ref_count_forward",
            f"{p}ref_count_reverse",
            f"{p}alt_count_forward",
            f"{p}alt_count_reverse",
            f"{p}ref_count_fragment_forward",
            f"{p}ref_count_fragment_reverse",
            f"{p}alt_count_fragment_forward",
            f"{p}alt_count_fragment_reverse",
        ]

    def _init_writer(self, original_headers: list[str]) -> None:
        """
        Initialize the CSV writer with dynamically constructed headers.

        Header order: original MAF columns first, then gbcms columns appended.
        Duplicate column names (already present in original) are skipped.

        Args:
            original_headers: Column names from the input MAF (or defaults for VCF→MAF).
        """
        gbcms_cols = self._gbcms_column_names()
        existing = set(original_headers)

        # Only append gbcms columns not already in the original headers
        new_cols = [c for c in gbcms_cols if c not in existing]
        self.fieldnames = list(original_headers) + new_cols

        self.writer = csv.DictWriter(
            self.file,
            fieldnames=self.fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        self.writer.writeheader()
        self._headers_written = True

        logger.debug(
            "MafWriter headers: %d original + %d gbcms = %d total columns",
            len(original_headers),
            len(new_cols),
            len(self.fieldnames),
        )

    def _populate_gbcms_counts(self, counts: Any) -> dict[str, str]:
        """
        Build the gbcms count columns dictionary with the configured prefix.

        Calculates VAF values and formats all count data as strings.

        Args:
            counts: BaseCounts object from the Rust engine.

        Returns:
            Dictionary mapping prefixed column names to string values.
        """
        p = self.column_prefix

        # Calculate VAFs with zero-division protection
        total_reads = counts.rd + counts.ad
        vaf = counts.ad / total_reads if total_reads > 0 else 0.0

        total_frags = counts.rdf + counts.adf
        vaf_frag = counts.adf / total_frags if total_frags > 0 else 0.0

        return {
            # Core counts
            f"{p}ref_count": str(counts.rd),
            f"{p}alt_count": str(counts.ad),
            f"{p}total_count": str(counts.dp),
            f"{p}vaf": f"{vaf:.4f}",
            # Fragment counts
            f"{p}ref_count_fragment": str(counts.rdf),
            f"{p}alt_count_fragment": str(counts.adf),
            f"{p}total_count_fragment": str(counts.dpf),
            f"{p}vaf_fragment": f"{vaf_frag:.4f}",
            # Strand bias (unprefixed)
            "strand_bias_p_value": f"{counts.sb_pval:.4e}",
            "strand_bias_odds_ratio": f"{counts.sb_or:.4f}",
            "fragment_strand_bias_p_value": f"{counts.fsb_pval:.4e}",
            "fragment_strand_bias_odds_ratio": f"{counts.fsb_or:.4f}",
            # Strand counts
            f"{p}ref_count_forward": str(counts.rd_fwd),
            f"{p}ref_count_reverse": str(counts.rd_rev),
            f"{p}alt_count_forward": str(counts.ad_fwd),
            f"{p}alt_count_reverse": str(counts.ad_rev),
            f"{p}ref_count_fragment_forward": str(counts.rdf_fwd),
            f"{p}ref_count_fragment_reverse": str(counts.rdf_rev),
            f"{p}alt_count_fragment_forward": str(counts.adf_fwd),
            f"{p}alt_count_fragment_reverse": str(counts.adf_rev),
        }

    def write(self, variant: Variant, counts: Any, sample_name: str = "TUMOR") -> None:
        """
        Write a single variant row to the MAF output.

        Two output strategies:
        - MAF→MAF (variant.metadata populated): Pass through all original columns,
          append gbcms count columns. Original values are NEVER overwritten.
        - VCF→MAF (no metadata): Generate GDC-compliant MAF coordinates from
          internal representation using CoordinateKernel.internal_to_maf().

        Args:
            variant: Normalized Variant with optional metadata from input MAF.
            counts: BaseCounts object from the Rust engine.
            sample_name: Sample name for Tumor_Sample_Barcode column.
        """
        # Initialize writer on first variant (headers depend on input format)
        if not self._headers_written:
            if variant.metadata:
                # MAF→MAF: use original input headers
                self._init_writer(list(variant.metadata.keys()))
            else:
                # VCF→MAF: use default GDC MAF headers + VCF-origin fields
                vcf_headers = self._DEFAULT_MAF_HEADERS + [
                    "vcf_id",
                    "vcf_pos",
                    "vcf_region",
                ]
                self._init_writer(vcf_headers)

        assert self.writer is not None

        # Build the output row based on input format
        if variant.metadata:
            # MAF→MAF: start with ALL original metadata (preserves every column)
            row = dict(variant.metadata)
        else:
            # VCF→MAF: build row from internal representation
            row = dict.fromkeys(self.fieldnames, "")

            # Convert internal coordinates to GDC MAF format
            maf_coords = CoordinateKernel.internal_to_maf(variant)
            row.update(maf_coords)
            row["Chromosome"] = variant.chrom

            # VCF-origin tracking fields
            vcf_pos = variant.pos + 1
            row["vcf_pos"] = str(vcf_pos)
            row["vcf_region"] = f"{variant.chrom}:{vcf_pos}"
            if variant.original_id:
                row["vcf_id"] = variant.original_id

        # Set sample barcode:
        # - VCF→MAF: always use BAM sample name (no barcode in VCF)
        # - MAF→MAF + preserve_barcode: keep original from input metadata
        # - MAF→MAF + no preserve_barcode: override with BAM sample name
        if not (variant.metadata and self.preserve_barcode):
            row["Tumor_Sample_Barcode"] = sample_name

        # Append gbcms count columns (both paths, never overwrites originals)
        row.update(self._populate_gbcms_counts(counts))

        self.writer.writerow(row)

    def close(self) -> None:
        """Close the output file."""
        self.file.close()
        logger.debug("MafWriter closed: %s", self.path)


class VcfWriter(OutputWriter):
    """Writes results to a VCF file."""

    def __init__(self, path: Path, sample_name: str = "SAMPLE"):
        self.path = path
        self.sample_name = sample_name
        self.file = open(path, "w")
        self._headers_written = False

    def _write_header(self):
        # Minimal VCF header
        headers = [
            "##fileformat=VCFv4.2",
            "##source=gbcms_v2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher strand bias p-value">',
            '##INFO=<ID=SB_OR,Number=1,Type=Float,Description="Fisher strand bias odds ratio">',
            '##INFO=<ID=FSB_PVAL,Number=1,Type=Float,Description="Fisher fragment strand bias p-value">',
            '##INFO=<ID=FSB_OR,Number=1,Type=Float,Description="Fisher fragment strand bias odds ratio">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles (fwd,rev)">',
            '##FORMAT=<ID=DP,Number=2,Type=Integer,Description="Approximate read depth (ref_total,alt_total)">',
            '##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference read depth (fwd,rev)">',
            '##FORMAT=<ID=RDF,Number=2,Type=Integer,Description="Ref Fragment Count (fwd,rev)">',
            '##FORMAT=<ID=ADF,Number=2,Type=Integer,Description="Alt Fragment Count (fwd,rev)">',
            '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction (read level)">',
            '##FORMAT=<ID=FAF,Number=1,Type=Float,Description="Variant Allele Fraction (fragment level)">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}",
        ]
        self.file.write("\n".join(headers) + "\n")
        self._headers_written = True

    def write(self, variant: Variant, counts: Any, sample_name: str = "SAMPLE"):
        if not self._headers_written:
            self._write_header()

        # VCF POS is 1-based
        pos = variant.pos + 1

        # INFO fields
        info = f"DP={counts.dp};SB_PVAL={counts.sb_pval:.4e};SB_OR={counts.sb_or:.4f};FSB_PVAL={counts.fsb_pval:.4e};FSB_OR={counts.fsb_or:.4f}"

        # FORMAT fields
        # GT: Simple 0/1 if alt > 0? Or ./1?
        # Let's assume 0/1 if we have alt counts, else 0/0
        gt = "0/1" if counts.ad > 0 else "0/0"

        # DP: ref_total,alt_total
        dp = f"{counts.rd},{counts.ad}"

        # RD: ref_fwd,ref_rev
        rd = f"{counts.rd_fwd},{counts.rd_rev}"

        # AD: alt_fwd,alt_rev
        ad = f"{counts.ad_fwd},{counts.ad_rev}"

        # RDF: ref_frag_fwd,ref_frag_rev
        rdf = f"{counts.rdf_fwd},{counts.rdf_rev}"

        # ADF: alt_frag_fwd,alt_frag_rev
        adf = f"{counts.adf_fwd},{counts.adf_rev}"

        # VAF calculations
        total_reads = counts.rd + counts.ad
        vaf = counts.ad / total_reads if total_reads > 0 else 0.0

        total_frags = counts.rdf + counts.adf
        faf = counts.adf / total_frags if total_frags > 0 else 0.0

        format_str = "GT:DP:RD:AD:RDF:ADF:VAF:FAF"
        sample_data = f"{gt}:{dp}:{rd}:{ad}:{rdf}:{adf}:{vaf:.4f}:{faf:.4f}"

        row = [
            variant.chrom,
            str(pos),
            variant.original_id or ".",
            variant.ref,
            variant.alt,
            ".",  # QUAL
            ".",  # FILTER
            info,
            format_str,
            sample_data,
        ]

        self.file.write("\t".join(row) + "\n")

    def close(self):
        self.file.close()
