"""
CLI Entry Point: Exposes the gbcms functionality via command line.
"""

import logging
from pathlib import Path

import typer

from .models.core import GbcmsConfig, OutputFormat
from .pipeline import Pipeline
from .utils import setup_logging

__all__ = ["app", "run"]

logger = logging.getLogger(__name__)

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
    output_suffix: str = typer.Option(
        "", "--suffix", "-S", help="Suffix to append to output filename (e.g. '.genotyped')"
    ),
    min_mapq: int = typer.Option(20, "--min-mapq", help="Minimum mapping quality"),
    min_baseq: int = typer.Option(0, "--min-baseq", help="Minimum base quality"),
    filter_duplicates: bool = typer.Option(True, help="Filter duplicate reads"),
    filter_secondary: bool = typer.Option(False, help="Filter secondary alignments"),
    filter_supplementary: bool = typer.Option(False, help="Filter supplementary alignments"),
    filter_qc_failed: bool = typer.Option(False, help="Filter reads failing QC"),
    filter_improper_pair: bool = typer.Option(False, help="Filter improperly paired reads"),
    filter_indel: bool = typer.Option(False, help="Filter reads containing indels"),
    threads: int = typer.Option(
        1, "--threads", "-t", help="Number of threads for parallel processing"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
):
    """
    Run gbcms on one or more BAM files.
    """
    # Configure logging
    setup_logging(verbose=verbose)

    # Parse BAM inputs
    bams_dict = _parse_bam_inputs(bam_files, bam_list)

    if not bams_dict:
        logger.error("No valid BAM files provided via --bam or --bam-list")
        raise typer.Exit(code=1)

    logger.info("Found %d BAM file(s) to process", len(bams_dict))

    try:
        config = GbcmsConfig(
            variant_file=variant_file,
            bam_files=bams_dict,
            reference_fasta=reference,
            output_dir=output_dir,
            output_format=output_format,
            output_suffix=output_suffix,
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
            filter_duplicates=filter_duplicates,
            filter_secondary=filter_secondary,
            filter_supplementary=filter_supplementary,
            filter_qc_failed=filter_qc_failed,
            filter_improper_pair=filter_improper_pair,
            filter_indel=filter_indel,
            threads=threads,
        )

        pipeline = Pipeline(config)
        pipeline.run()

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)
        raise typer.Exit(code=1) from e


def _parse_bam_inputs(
    bam_files: list[Path] | None, bam_list: Path | None
) -> dict[str, Path]:
    """
    Parse BAM inputs from direct arguments and/or BAM list file.

    Args:
        bam_files: List of BAM paths (optionally with sample_id:path format).
        bam_list: Path to file containing BAM paths (one per line).

    Returns:
        Dictionary mapping sample names to BAM paths.
    """
    bams_dict: dict[str, Path] = {}

    # 1. Process direct BAM arguments
    if bam_files:
        for bam_arg in bam_files:
            sample_name, bam_path = _parse_bam_arg(bam_arg)

            if not bam_path.exists():
                logger.error("BAM file not found: %s", bam_path)
                continue

            bams_dict[sample_name] = bam_path

    # 2. Process BAM list file
    if bam_list:
        if not bam_list.exists():
            logger.error("BAM list file not found: %s", bam_list)
            return bams_dict

        try:
            with open(bam_list) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split()
                    if len(parts) >= 2:
                        sample_name = parts[0]
                        bam_path = Path(parts[1])
                    else:
                        bam_path = Path(parts[0])
                        sample_name = bam_path.stem

                    if not bam_path.exists():
                        logger.warning("BAM file from list not found: %s", bam_path)
                        continue

                    bams_dict[sample_name] = bam_path

        except Exception as e:
            logger.error("Error reading BAM list file %s: %s", bam_list, e)

    return bams_dict


def _parse_bam_arg(bam_arg: Path) -> tuple[str, Path]:
    """
    Parse a BAM argument that may be in sample_id:path format.

    Args:
        bam_arg: Path object (may contain sample_id:path as string).

    Returns:
        Tuple of (sample_name, bam_path).
    """
    bam_str = str(bam_arg)
    if ":" in bam_str:
        parts = bam_str.split(":", 1)
        return parts[0], Path(parts[1])
    return bam_arg.stem, bam_arg


if __name__ == "__main__":
    app()
