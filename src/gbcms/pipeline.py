"""
Pipeline Orchestrator: Manages the execution flow of gbcms.

This module handles:
1. Reading variants from input (VCF/MAF).
2. Iterating over samples (BAM files).
3. Running the Rust-based counting engine for each sample.
4. Writing results to per-sample output files.
"""

from pathlib import Path
from typing import List, Dict, Optional, Union
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


class Pipeline:
    def __init__(self, config: GbcmsConfig):
        self.config = config
        self.console = Console()

    def run(self):
        """Execute the pipeline."""
        self.console.print("[bold blue]Starting gbcms pipeline[/bold blue]")
        self.console.print(f"Output directory: {self.config.output_dir}")

        # 1. Load Variants
        with self.console.status("[bold green]Loading variants...[/bold green]"):
            variants = self._load_variants()

        self.console.print(f"Loaded [bold]{len(variants)}[/bold] variants.")

        if not variants:
            self.console.print("[bold red]No variants found. Exiting.[/bold red]")
            return

        # 2. Validate Variants against Reference
        with self.console.status(
            "[bold green]Validating variants against reference...[/bold green]"
        ):
            valid_variants = self._validate_variants(variants)

        self.console.print(f"Valid variants: [bold]{len(valid_variants)}[/bold] / {len(variants)}")

        if not valid_variants:
            self.console.print(
                "[bold red]No valid variants remaining after validation. Exiting.[/bold red]"
            )
            return

        variants = valid_variants

        # 3. Prepare Rust Variants
        rs_variants = [
            gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants
        ]

        # 4. Process Each Sample
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

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

                # Validate BAM Header
                if not self._validate_bam_header(bam_path, variants):
                    self.console.print(
                        f"[yellow]Warning: BAM {sample_name} may not contain variant chromosomes. Proceeding anyway...[/yellow]"
                    )

                try:
                    # Run Rust Engine
                    counts_list = gbcms_rs.count_bam(
                        str(bam_path),
                        rs_variants,
                        min_mapq=self.config.min_mapping_quality,
                        min_baseq=self.config.min_base_quality,
                        filter_duplicates=self.config.filter_duplicates,
                        filter_secondary=self.config.filter_secondary,
                        filter_supplementary=self.config.filter_supplementary,
                    )

                    # Write Output
                    self._write_output(sample_name, variants, counts_list)

                except Exception as e:
                    self.console.print(
                        f"[bold red]Error processing sample {sample_name}: {e}[/bold red]"
                    )
                    # Continue to next sample

                progress.advance(task)

        self.console.print("[bold green]Pipeline completed successfully.[/bold green]")

    def _load_variants(self) -> list[Variant]:
        """Load variants based on file extension."""
        path = self.config.variant_file
        reader: VariantReader

        if path.suffix.lower() in [".vcf", ".gz"]:  # .vcf.gz handled by pysam
            reader = VcfReader(path)
        elif path.suffix.lower() == ".maf":
            reader = MafReader(path, fasta_path=self.config.reference_fasta)
        else:
            raise ValueError(f"Unsupported variant file format: {path.suffix}")

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
                if invalid_count <= 5:  # Log first few failures
                    self.console.print(
                        f"[yellow]Invalid variant (REF mismatch): {v.chrom}:{v.pos} {v.ref}>{v.alt}[/yellow]"
                    )

        if invalid_count > 5:
            self.console.print(
                f"[yellow]... and {invalid_count - 5} more invalid variants.[/yellow]"
            )

        checker.close()
        return valid_variants

    def _validate_bam_header(self, bam_path: Path, variants: list[Variant]) -> bool:
        """Check if BAM header contains chromosomes from variants."""
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                bam_chroms = set(bam.references)

            # Check a few variants
            # We need to handle chr prefix normalization
            # BAM might have 'chr1', variant '1', or vice versa.

            # Normalize BAM chroms
            norm_bam_chroms = {CoordinateKernel.normalize_chromosome(c) for c in bam_chroms}

            # Check first variant as a heuristic
            if variants:
                v = variants[0]
                norm_v_chrom = CoordinateKernel.normalize_chromosome(v.chrom)
                if norm_v_chrom not in norm_bam_chroms:
                    return False
            return True
        except Exception as e:
            self.console.print(f"[yellow]Could not validate BAM header: {e}[/yellow]")
            return True  # Assume ok if we can't check

    def _write_output(
        self, sample_name: str, variants: list[Variant], counts_list: list[gbcms_rs.BaseCounts]
    ):
        """Write results to output file."""
        ext = "vcf" if self.config.output_format == OutputFormat.VCF else "maf"
        output_path = self.config.output_dir / f"{sample_name}.{ext}" # Re-added original path construction
        writer: Union[VcfWriter, MafWriter]
        if self.config.output_format == OutputFormat.VCF:
            writer = VcfWriter(output_path, sample_name=sample_name)
        else:
            writer = MafWriter(output_path)
        
        for v, counts in zip(variants, counts_list): # Changed 'results' to 'counts_list'
            writer.write(v, counts, sample_name=sample_name)
            
        writer.close()

        # self.console.print(f"Results written to {output_path}")
