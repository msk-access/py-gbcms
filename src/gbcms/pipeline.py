"""
Pipeline Orchestrator: Manages the execution flow of gbcms.

This module handles:
1. Reading variants from input (VCF/MAF).
2. Iterating over samples (BAM files).
3. Running the Rust-based counting engine for each sample.
4. Writing results to per-sample output files.
"""

import logging
import time
from pathlib import Path

import pysam
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
)

import gbcms_rs

from .core.kernel import CoordinateKernel
from .io.input import MafReader, ReferenceChecker, VariantReader, VcfReader
from .io.output import MafWriter, VcfWriter
from .models.core import GbcmsConfig, OutputFormat, Variant

logger = logging.getLogger(__name__)

__all__ = ["Pipeline"]


class Pipeline:
    """Main pipeline for processing BAM files and counting bases at variant positions."""

    def __init__(self, config: GbcmsConfig):
        """
        Initialize the pipeline.

        Args:
            config: Configuration object with input/output paths and filter settings.
        """
        self.config = config
        self.console = Console()
        self._stats = {"samples_processed": 0, "total_variants": 0, "total_time": 0.0}

    def run(self) -> dict:
        """
        Execute the pipeline.

        Returns:
            Dictionary with processing statistics.
        """
        start_time = time.perf_counter()
        logger.info("Starting gbcms pipeline")
        logger.info("Output directory: %s", self.config.output.directory)

        # 1. Load Variants
        logger.debug("Loading variants from %s", self.config.variant_file)
        variants = self._load_variants()
        logger.info("Loaded %d variants", len(variants))

        if not variants:
            logger.error("No variants found. Exiting.")
            return self._stats

        # 2. Validate Variants against Reference
        logger.debug("Validating variants against reference genome")
        valid_variants = self._validate_variants(variants)
        logger.info("Valid variants: %d / %d", len(valid_variants), len(variants))

        if not valid_variants:
            logger.error("No valid variants remaining after validation. Exiting.")
            return self._stats

        variants = valid_variants
        self._stats["total_variants"] = len(variants)

        # 3. Prepare Rust Variants
        rs_variants = [
            gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value)
            for v in variants
        ]

        # 4. Process Each Sample
        self.config.output.directory.mkdir(parents=True, exist_ok=True)
        samples = list(self.config.bam_files.items())

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
            console=self.console,
        ) as progress:
            task = progress.add_task("[cyan]Processing samples...", total=len(samples))

            for sample_name, bam_path in samples:
                progress.update(task, description=f"[cyan]Processing {sample_name}...")
                self._process_sample(sample_name, bam_path, variants, rs_variants)
                progress.advance(task)

        # Calculate total time
        self._stats["total_time"] = time.perf_counter() - start_time
        logger.info(
            "Pipeline completed: %d samples, %.2fs",
            self._stats["samples_processed"],
            self._stats["total_time"],
        )

        return self._stats

    def _process_sample(
        self,
        sample_name: str,
        bam_path: Path,
        variants: list[Variant],
        rs_variants: list,
    ) -> None:
        """
        Process a single sample.

        Args:
            sample_name: Name of the sample.
            bam_path: Path to BAM file.
            variants: List of normalized variants.
            rs_variants: List of Rust variant objects.
        """
        sample_start = time.perf_counter()
        logger.debug("Processing sample: %s (%s)", sample_name, bam_path)

        # Validate BAM Header
        if not self._validate_bam_header(bam_path, variants):
            logger.warning(
                "BAM %s may not contain variant chromosomes. Proceeding anyway.",
                sample_name,
            )

        try:
            # Run Rust Engine with nested config accessors
            rust_start = time.perf_counter()
            counts_list = gbcms_rs.count_bam(
                str(bam_path),
                rs_variants,
                min_mapq=self.config.quality.min_mapping_quality,
                min_baseq=self.config.quality.min_base_quality,
                filter_duplicates=self.config.filters.duplicates,
                filter_secondary=self.config.filters.secondary,
                filter_supplementary=self.config.filters.supplementary,
                filter_qc_failed=self.config.filters.qc_failed,
                filter_improper_pair=self.config.filters.improper_pair,
                filter_indel=self.config.filters.indel,
                threads=self.config.threads,
            )
            rust_time = time.perf_counter() - rust_start
            logger.debug("Rust count_bam completed in %.3fs", rust_time)

            # Write Output
            self._write_output(sample_name, variants, counts_list)
            self._stats["samples_processed"] += 1

            sample_time = time.perf_counter() - sample_start
            logger.debug("Sample %s completed in %.3fs", sample_name, sample_time)

        except Exception as e:
            logger.error("Error processing sample %s: %s", sample_name, e)

    def _load_variants(self) -> list[Variant]:
        """Load variants based on file extension."""
        path = self.config.variant_file
        reader: VariantReader

        suffix = path.suffix.lower()
        if suffix in [".vcf", ".gz"]:
            reader = VcfReader(path)
        elif suffix == ".maf":
            reader = MafReader(path, fasta_path=self.config.reference_fasta)
        else:
            raise ValueError(f"Unsupported variant file format: {suffix}")

        variants = list(reader)
        if hasattr(reader, "close"):
            reader.close()

        return variants

    def _validate_variants(self, variants: list[Variant]) -> list[Variant]:
        """Validate variants against reference genome."""
        checker = ReferenceChecker(self.config.reference_fasta)
        valid_variants = []
        invalid_count = 0

        for v in variants:
            if checker.validate(v):
                valid_variants.append(v)
            else:
                invalid_count += 1
                if invalid_count <= 5:
                    logger.warning(
                        "Invalid variant (REF mismatch): %s:%d %s>%s",
                        v.chrom,
                        v.pos,
                        v.ref,
                        v.alt,
                    )

        if invalid_count > 5:
            logger.warning("... and %d more invalid variants.", invalid_count - 5)

        checker.close()
        return valid_variants

    def _validate_bam_header(self, bam_path: Path, variants: list[Variant]) -> bool:
        """Check if BAM header contains chromosomes from variants."""
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                bam_chroms = set(bam.references)

            norm_bam_chroms = {
                CoordinateKernel.normalize_chromosome(c) for c in bam_chroms
            }

            if variants:
                v = variants[0]
                norm_v_chrom = CoordinateKernel.normalize_chromosome(v.chrom)
                if norm_v_chrom not in norm_bam_chroms:
                    return False
            return True
        except Exception as e:
            logger.warning("Could not validate BAM header: %s", e)
            return True

    def _write_output(
        self,
        sample_name: str,
        variants: list[Variant],
        counts_list: list[gbcms_rs.BaseCounts],
    ) -> None:
        """Write results to output file."""
        ext = "vcf" if self.config.output.format == OutputFormat.VCF else "maf"
        suffix = self.config.output.suffix
        output_path = self.config.output.directory / f"{sample_name}{suffix}.{ext}"

        writer: VcfWriter | MafWriter
        if self.config.output.format == OutputFormat.VCF:
            writer = VcfWriter(output_path, sample_name=sample_name)
        else:
            writer = MafWriter(output_path)

        for v, counts in zip(variants, counts_list, strict=True):
            writer.write(v, counts, sample_name=sample_name)

        writer.close()
        logger.debug("Results written to %s", output_path)
