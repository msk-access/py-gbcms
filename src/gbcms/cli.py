"""
CLI Entry Point: Exposes the gbcms functionality via command line.
"""

from pathlib import Path

import typer

from .models.core import GbcmsConfig, OutputFormat
from .pipeline import Pipeline

app = typer.Typer(help="gbcms: Get Base Counts Multi-Sample")


@app.callback()
def main():
    """
    gbcms: Get Base Counts Multi-Sample
    """
    pass


@app.command()
def run(
    variant_file: Path = typer.Option(
        ..., "--variants", "-v", help="Path to VCF or MAF file containing variants"
    ),
    bam_files: list[Path] | None = typer.Option(
        None, "--bam", "-b", help="Path to BAM file(s). Can be specified multiple times."
    ),
    bam_list: Path | None = typer.Option(
        None, "--bam-list", "-L", help="File containing list of BAM paths (one per line)"
    ),
    reference: Path = typer.Option(..., "--fasta", "-f", help="Path to reference FASTA file"),
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o", help="Directory to write output files"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.VCF, "--format", help="Output format (vcf or maf)"
    ),
    min_mapq: int = typer.Option(20, "--min-mapq", help="Minimum mapping quality"),
    min_baseq: int = typer.Option(0, "--min-baseq", help="Minimum base quality"),
    filter_duplicates: bool = typer.Option(True, help="Filter duplicate reads"),
    threads: int = typer.Option(
        1, "--threads", "-t", help="Number of threads (not yet implemented in v2 python layer)"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-V", help="Enable verbose debug logging"
    ),
):
    """
    Run gbcms on one or more BAM files.
    """
    from rich.console import Console
    from rich.logging import RichHandler
    import logging

    # Configure logging
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, markup=True)]
    )

    console = Console()

    # Map BAMs to sample names (filename stem for now)
    bams_dict = {}

    # 1. Process direct BAM arguments
    if bam_files:
        for bam_path in bam_files:
            if not bam_path.exists():
                console.print(f"[bold red]Error: BAM file not found: {bam_path}[/bold red]")
                raise typer.Exit(code=1)
            sample_name = bam_path.stem
            bams_dict[sample_name] = bam_path

    # 2. Process BAM list file
    if bam_list:
        if not bam_list.exists():
            console.print(f"[bold red]Error: BAM list file not found: {bam_list}[/bold red]")
            raise typer.Exit(code=1)

        try:
            with open(bam_list) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    bam_path = Path(line)
                    if not bam_path.exists():
                        console.print(
                            f"[yellow]Warning: BAM file from list not found: {bam_path}[/yellow]"
                        )
                        continue
                    sample_name = bam_path.stem
                    bams_dict[sample_name] = bam_path
        except Exception as e:
            console.print(f"[bold red]Error reading BAM list file {bam_list}: {e}[/bold red]")
            raise typer.Exit(code=1) from e

    if not bams_dict:
        console.print(
            "[bold red]Error: No valid BAM files provided via --bam or --bam-list[/bold red]"
        )
        raise typer.Exit(code=1)

    try:
        config = GbcmsConfig(
            variant_file=variant_file,
            bam_files=bams_dict,
            reference_fasta=reference,
            output_dir=output_dir,
            output_format=output_format,
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
            filter_duplicates=filter_duplicates,
            threads=threads,
        )

        pipeline = Pipeline(config)
        pipeline.run()

    except Exception as e:
        console.print(f"[bold red]Error: {e}[/bold red]")
        raise typer.Exit(code=1) from e


if __name__ == "__main__":
    app()
