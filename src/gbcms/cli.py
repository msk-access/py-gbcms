"""
CLI Entry Point: Exposes the gbcms functionality via command line.

Validation order (enforced for every option):
  1. Parse-time  — Typer Enum / min=/max= for constrained choices and numeric ranges.
  2. Pre-model   — Explicit checks in the command body before Pydantic construction
                   (file extensions, cross-option semantics, charset validation).
  3. Model-time  — Pydantic field constraints and validators in models/core.py.
  4. No silent skips — Missing inputs fail-fast unless the caller opts-out explicitly
                       (e.g. --lenient-bam).
"""

import logging
import os
import re
from pathlib import Path

import typer

from . import __version__
from .models.core import (
    AlignmentConfig,
    GbcmsConfig,
    OutputConfig,
    OutputFormat,
    QualityThresholds,
    ReadFilters,
    StrEnum,  # canonical backport (Python ≮3.10 compatible), defined in models.core
)
from .pipeline import Pipeline
from .utils import setup_logging

__all__ = ["app", "run", "normalize"]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constrained-choice enums (parse-time validation via Typer)
# ---------------------------------------------------------------------------


class AlignmentBackend(StrEnum):
    """CLI-exposed alignment backend options.

    'pairhmm' is accepted as a model-level alias but intentionally not
    surfaced here so users always type the canonical short form.
    """

    SW = "sw"
    HMM = "hmm"


# Valid variant file extensions (checked before Pydantic config construction)
_VALID_VARIANT_EXTENSIONS: frozenset[str] = frozenset({".vcf", ".maf"})

# Column-prefix charset: only letters, digits, underscores
_COLUMN_PREFIX_RE = re.compile(r"^[A-Za-z0-9_]*$")

app = typer.Typer(help="gbcms: Get Base Counts Multi-Sample")


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        typer.echo(f"gbcms {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool | None = typer.Option(
        None, "--version", callback=version_callback, is_eager=True, help="Show version and exit."
    ),
) -> None:
    """
    gbcms: Get Base Counts Multi-Sample
    """
    pass


@app.command()
def run(
    # Input options
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
    # Output options
    output_dir: Path = typer.Option(
        ..., "--output-dir", "-o", help="Directory to write output files"
    ),
    output_format: OutputFormat = typer.Option(
        OutputFormat.VCF, "--format", help="Output format (vcf or maf)"
    ),
    output_suffix: str = typer.Option(
        "", "--suffix", "-S", help="Suffix to append to output filename (e.g. '.genotyped')"
    ),
    column_prefix: str = typer.Option(
        "",
        "--column-prefix",
        help=(
            "Prefix for gbcms count columns in MAF output. "
            "Default: no prefix (e.g., 'ref_count'). "
            "Use 't_' for legacy compatibility (e.g., 't_ref_count')."
        ),
    ),
    preserve_barcode: bool = typer.Option(
        False,
        "--preserve-barcode",
        help=(
            "Preserve original Tumor_Sample_Barcode from input MAF "
            "instead of overriding with BAM sample name. "
            "Only applies to MAF→MAF output."
        ),
    ),
    # Quality thresholds
    min_mapq: int = typer.Option(20, "--min-mapq", help="Minimum mapping quality"),
    min_baseq: int = typer.Option(20, "--min-baseq", help="Minimum base quality"),
    fragment_qual_threshold: int = typer.Option(
        10,
        "--fragment-qual-threshold",
        help=(
            "Quality difference threshold for fragment consensus. "
            "When R1 and R2 disagree, the higher-quality allele wins only if "
            "the difference exceeds this threshold; otherwise the fragment is discarded."
        ),
    ),
    context_padding: int = typer.Option(
        5,
        "--context-padding",
        min=1,
        max=50,
        help=(
            "Minimum flanking reference bases around indel/complex variants for "
            "haplotype construction and SW alignment. Range 1–50 enforced at "
            "parse time. Auto-increased in repeat regions when --adaptive-context is enabled."
        ),
    ),
    adaptive_context: bool = typer.Option(
        True,
        "--adaptive-context/--no-adaptive-context",
        help="Dynamically increase context padding in tandem repeat regions.",
    ),
    # Read filters
    filter_duplicates: bool = typer.Option(True, help="Filter duplicate reads"),
    filter_secondary: bool = typer.Option(False, help="Filter secondary alignments"),
    filter_supplementary: bool = typer.Option(False, help="Filter supplementary alignments"),
    filter_qc_failed: bool = typer.Option(False, help="Filter reads failing QC"),
    filter_improper_pair: bool = typer.Option(False, help="Filter improperly paired reads"),
    filter_indel: bool = typer.Option(False, help="Filter reads containing indels"),
    # Normalization
    show_normalization: bool = typer.Option(
        False,
        "--show-normalization",
        help="Add norm_* columns showing left-aligned coordinates in output.",
    ),
    # Performance
    threads: int = typer.Option(
        1, "--threads", "-t", help="Number of threads for parallel processing"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
    trace: bool = typer.Option(
        False,
        "--trace",
        "-T",
        help="Enable per-read Rust trace logging (slow). Implies --verbose. "
        "Shows detailed per-read classification diagnostics from the counting engine.",
    ),
    # BAM robustness
    lenient_bam: bool = typer.Option(
        False,
        "--lenient-bam",
        help=(
            "Skip missing BAM files and continue with the remaining samples. "
            "Default (off): any missing BAM causes an immediate exit. "
            "Use this flag when running on batch lists where some files may be absent."
        ),
    ),
    # Alignment backend (advanced)
    alignment_backend: AlignmentBackend = typer.Option(
        AlignmentBackend.SW,
        "--alignment-backend",
        help=(
            "Alignment backend for Phase 3 classification: "
            "'sw' (Smith-Waterman, default) or 'hmm' (PairHMM). "
            "Invalid values are rejected at parse time."
        ),
    ),
    hmm_llr_threshold: float = typer.Option(
        2.3,
        "--llr-threshold",
        help="PairHMM log-likelihood ratio threshold for confident calls (default: ln(10) ≈ 2.3).",
    ),
    hmm_gap_open: float = typer.Option(
        1e-4,
        "--gap-open-prob",
        help="PairHMM gap-open probability for non-repeat regions.",
    ),
    hmm_gap_extend: float = typer.Option(
        0.1,
        "--gap-extend-prob",
        help="PairHMM gap-extend probability for non-repeat regions.",
    ),
    hmm_gap_open_repeat: float = typer.Option(
        1e-2,
        "--repeat-gap-open-prob",
        help="PairHMM gap-open probability for tandem repeat regions.",
    ),
    hmm_gap_extend_repeat: float = typer.Option(
        0.5,
        "--repeat-gap-extend-prob",
        help="PairHMM gap-extend probability for tandem repeat regions.",
    ),
):
    """
    Run gbcms on one or more BAM files.
    """
    # ── 1. Logging (must be first so all subsequent checks log correctly) ──────
    setup_logging(verbose=verbose, trace=trace)

    # ── 2. Pre-model validation (semantic + cross-option checks) ──────────────

    # GAP 12: Reject unsupported variant file extensions before any I/O.
    _is_vcf_gz = variant_file.name.lower().endswith(".vcf.gz")
    _ext = variant_file.suffix.lower()
    if not _is_vcf_gz and _ext not in _VALID_VARIANT_EXTENSIONS:
        logger.error(
            "Unsupported variant file extension '%s'. Expected .vcf, .vcf.gz, or .maf. "
            "Got: %s",
            _ext,
            variant_file,
        )
        raise typer.Exit(code=1)

    # GAP 10: Validate --column-prefix charset (letters, digits, underscores only).
    if column_prefix and not _COLUMN_PREFIX_RE.match(column_prefix):
        logger.error(
            "Invalid --column-prefix '%s': only letters, digits, and underscores are allowed. "
            "Whitespace and special characters would produce malformed column names.",
            column_prefix,
        )
        raise typer.Exit(code=1)

    # GAP 9: Warn when --preserve-barcode is used with non-MAF input (it is a no-op).
    if preserve_barcode and not _is_vcf_gz and _ext != ".maf":
        logger.warning(
            "--preserve-barcode has no effect when the variant file is not a MAF "
            "(got '%s'). The BAM sample name will be used in all output rows.",
            variant_file.suffix,
        )

    # GAP 8: Advisory warning when threads exceeds available CPUs.
    cpu_count = os.cpu_count() or 1
    if threads > cpu_count:
        logger.warning(
            "--threads %d exceeds os.cpu_count() (%d). "
            "Performance may degrade due to CPU oversubscription.",
            threads,
            cpu_count,
        )

    # ── 3. Parse BAM inputs (fail-fast by default; --lenient-bam opts out) ─────
    bams_dict = _parse_bam_inputs(bam_files, bam_list, lenient=lenient_bam)

    if not bams_dict:
        logger.error(
            "No BAM files to process. Provide at least one via --bam <path> "
            "or --bam-list <file>. Use --lenient-bam to allow partial BAM lists."
        )
        raise typer.Exit(code=1)

    logger.info("Found %d BAM file(s) to process", len(bams_dict))

    try:
        # Build nested config objects
        output_config = OutputConfig(
            directory=output_dir,
            format=output_format,
            suffix=output_suffix,
            column_prefix=column_prefix,
            preserve_barcode=preserve_barcode,
        )

        quality_config = QualityThresholds(
            min_mapping_quality=min_mapq,
            min_base_quality=min_baseq,
            fragment_qual_threshold=fragment_qual_threshold,
            context_padding=context_padding,
            adaptive_context=adaptive_context,
        )

        filter_config = ReadFilters(
            duplicates=filter_duplicates,
            secondary=filter_secondary,
            supplementary=filter_supplementary,
            qc_failed=filter_qc_failed,
            improper_pair=filter_improper_pair,
            indel=filter_indel,
        )

        # Pass .value so AlignmentConfig receives a plain str, not the enum wrapper.
        # This is required because AlignmentConfig.validate_backend operates on str.
        alignment_config = AlignmentConfig(
            backend=alignment_backend.value,
            hmm_llr_threshold=hmm_llr_threshold,
            hmm_gap_open=hmm_gap_open,
            hmm_gap_extend=hmm_gap_extend,
            hmm_gap_open_repeat=hmm_gap_open_repeat,
            hmm_gap_extend_repeat=hmm_gap_extend_repeat,
        )

        config = GbcmsConfig(
            variant_file=variant_file,
            bam_files=bams_dict,
            reference_fasta=reference,
            output=output_config,
            quality=quality_config,
            filters=filter_config,
            threads=threads,
            alignment=alignment_config,
            show_normalization=show_normalization,
        )

        pipeline = Pipeline(config)
        pipeline.run()

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)
        raise typer.Exit(code=1) from e


@app.command()
def normalize(
    variant_file: Path = typer.Option(
        ..., "--variants", "-v", help="Path to VCF or MAF file containing variants"
    ),
    reference: Path = typer.Option(..., "--fasta", "-f", help="Path to reference FASTA file"),
    output: Path = typer.Option(
        ..., "--output", "-o", help="Output file path (TSV with normalization results)"
    ),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
    verbose: bool = typer.Option(False, "--verbose", "-V", help="Enable verbose debug logging"),
    trace: bool = typer.Option(
        False,
        "--trace",
        "-T",
        help="Enable per-read Rust trace logging (slow). Implies --verbose.",
    ),
):
    """
    Normalize variants (left-align + validate REF) without counting.

    Reads variants from a VCF or MAF file, applies MAF anchor resolution,
    REF validation, and bcftools-style left-alignment, then writes results
    to a TSV file showing original and normalized coordinates.
    """
    from .normalize import normalize_variants

    setup_logging(verbose=verbose, trace=trace)

    # Apply the same file extension pre-check as the 'run' command.
    _is_vcf_gz = variant_file.name.lower().endswith(".vcf.gz")
    _ext = variant_file.suffix.lower()
    if not _is_vcf_gz and _ext not in _VALID_VARIANT_EXTENSIONS:
        logger.error(
            "Unsupported variant file extension '%s'. Expected .vcf, .vcf.gz, or .maf. "
            "Got: %s",
            _ext,
            variant_file,
        )
        raise typer.Exit(code=1)

    normalize_variants(
        variant_file=variant_file,
        reference=reference,
        output=output,
        threads=threads,
    )


def _parse_bam_inputs(
    bam_files: list[Path] | None,
    bam_list: Path | None,
    *,
    lenient: bool = False,
) -> dict[str, Path]:
    """
    Parse BAM inputs from direct arguments and/or a BAM list file.

    Validation behaviour:
    - **Fail-fast (default)**: If any BAM path does not exist, all missing paths
      are logged at ERROR level and ``typer.Exit(code=1)`` is raised.
    - **Lenient mode** (``lenient=True``, enabled via ``--lenient-bam``): Missing
      paths are logged as errors but skipped; the run continues with the
      remaining samples.
    - **BAM list file not found**: Always fails immediately regardless of lenient
      mode.  The list file itself is a required input, not an optional sample.

    Args:
        bam_files: List of BAM paths (optionally with ``sample_id:path`` format).
        bam_list: Path to a file containing BAM paths (one per line,
            optionally ``sample_name<whitespace>path``).
        lenient: When True, skip missing BAM files instead of exiting.

    Returns:
        Dictionary mapping sample names to resolved BAM ``Path`` objects.

    Raises:
        typer.Exit: If any BAM file or the list file itself is missing and
            ``lenient`` is False.
    """
    bams_dict: dict[str, Path] = {}

    # ── 1. Process direct --bam arguments ────────────────────────────────────
    if bam_files:
        missing: list[str] = []
        for bam_arg in bam_files:
            sample_name, bam_path = _parse_bam_arg(bam_arg)

            if not bam_path.exists():
                logger.error("BAM file not found: %s", bam_path)
                missing.append(str(bam_path))
                continue

            logger.debug("Registered BAM sample '%s': %s", sample_name, bam_path)
            bams_dict[sample_name] = bam_path

        if missing and not lenient:
            logger.error(
                "%d BAM file(s) not found. "
                "Add --lenient-bam to skip missing files and continue with the rest.",
                len(missing),
            )
            raise typer.Exit(code=1)

    # ── 2. Process --bam-list file ────────────────────────────────────────────
    if bam_list:
        # The list file itself is always required — lenient mode does not apply here.
        if not bam_list.exists():
            logger.error(
                "BAM list file not found: %s. "
                "Note: --lenient-bam does not apply to the list file itself.",
                bam_list,
            )
            raise typer.Exit(code=1)

        logger.debug("Reading BAM list from: %s", bam_list)
        try:
            with open(bam_list) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue  # skip blanks and comment lines

                    parts = line.split()
                    if len(parts) >= 2:
                        sample_name = parts[0]
                        bam_path = Path(parts[1])
                    else:
                        bam_path = Path(parts[0])
                        sample_name = bam_path.stem

                    if not bam_path.exists():
                        # Upgraded from WARNING to ERROR — a missing BAM in the
                        # list is always unexpected, whether in lenient mode or not.
                        logger.error(
                            "BAM file from list not found: %s (sample '%s')",
                            bam_path,
                            sample_name,
                        )
                        continue

                    logger.debug(
                        "Registered BAM sample from list '%s': %s", sample_name, bam_path
                    )
                    bams_dict[sample_name] = bam_path

        except OSError as e:
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
