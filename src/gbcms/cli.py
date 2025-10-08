from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from . import __version__
from .config import Config
from .processor import VariantProcessor

# Initialize Typer app with rich help
app = typer.Typer(
    name="gbcms",
    help="Python implementation of gbcms for calculating base counts in BAM files",
    add_completion=False,
    rich_markup_mode="rich",
    no_args_is_help=True,
)

# Initialize Rich console
console = Console()


# Subcommands
count_app = typer.Typer(help="Count bases at variant positions")
validate_app = typer.Typer(help="Validate input files")
app.add_typer(count_app, name="count")
app.add_typer(validate_app, name="validate")


def validate_input_files(
    fasta: Path,
    bam_files: dict[str, str],
    variant_files: list[str],
    input_is_maf: bool,
    input_is_vcf: bool,
    rich_output: bool = False,
) -> tuple[bool, Table | None]:
    """
    Validate input files for gbcms processing.

    Args:
        fasta: Path to reference FASTA file
        bam_files: Dictionary of sample names to BAM file paths
        variant_files: List of variant file paths
        input_is_maf: Whether input files are in MAF format
        input_is_vcf: Whether input files are in VCF format
        rich_output: Whether to return detailed rich table for visual output

    Returns:
        Tuple of (is_valid, results_table_or_None)
    """
    from pathlib import Path

    if rich_output:
        from rich.table import Table

        results = Table(title="Validation Results", show_header=True, header_style="bold cyan")
        results.add_column("File Type", style="cyan")
        results.add_column("File Path", style="white")
        results.add_column("Status", style="white")
        results.add_column("Details", style="yellow")

        console.print(
            Panel.fit(
                "[bold cyan]File Validation[/bold cyan]\n" "Checking input files for gbcms",
                border_style="cyan",
            )
        )

    all_valid = True

    # Validate FASTA file and index
    if not fasta or str(fasta) == "":
        # Skip FASTA validation if not provided
        pass
    elif not fasta.exists():
        if rich_output:
            results.add_row("FASTA", str(fasta), "❌ FAIL", "File not found")
        else:
            console.print(f"[red]Error:[/red] FASTA file not found: {fasta}")
        all_valid = False
    else:
        fai_file = Path(str(fasta) + ".fai")
        if not fai_file.exists():
            if rich_output:
                results.add_row("FASTA", str(fasta), "⚠️  WARN", "Index (.fai) not found")
            else:
                console.print(f"[red]Error:[/red] FASTA index (.fai) not found: {fai_file}")
                console.print("Please index your FASTA file with: samtools faidx reference.fa")
            all_valid = False
        else:
            if rich_output:
                results.add_row("FASTA", str(fasta), "✅ PASS", "File and index found")

    # Validate BAM files and indices
    for sample_name, bam_path in bam_files.items():
        bam_file = Path(bam_path)
        if not bam_file.exists():
            if rich_output:
                results.add_row("BAM", f"{sample_name}:{bam_path}", "❌ FAIL", "File not found")
            else:
                console.print(f"[red]Error:[/red] BAM file not found: {bam_file}")
            all_valid = False
        else:
            # Check for BAM index files
            bai_file1 = Path(str(bam_file).replace(".bam", ".bai"))
            bai_file2 = Path(str(bam_file) + ".bai")

            if not bai_file1.exists() and not bai_file2.exists():
                if rich_output:
                    results.add_row(
                        "BAM", f"{sample_name}:{bam_path}", "⚠️  WARN", "Index (.bai) not found"
                    )
                else:
                    console.print(f"[red]Error:[/red] BAM index not found for: {bam_file}")
                    console.print(f"Expected: {bai_file1} or {bai_file2}")
                    console.print("Please index your BAM file with: samtools index sample.bam")
                all_valid = False
            else:
                if rich_output:
                    results.add_row(
                        "BAM", f"{sample_name}:{bam_path}", "✅ PASS", "File and index found"
                    )

    # Validate variant files
    for variant_file in variant_files:
        vcf = Path(variant_file)
        if not vcf.exists():
            if rich_output:
                results.add_row(
                    "VCF" if input_is_vcf else "MAF", str(vcf), "❌ FAIL", "File not found"
                )
            else:
                console.print(f"[red]Error:[/red] Variant file not found: {vcf}")
            all_valid = False
        else:
            if rich_output:
                results.add_row("VCF" if input_is_vcf else "MAF", str(vcf), "✅ PASS", "File found")

    if rich_output:
        return all_valid, results
    else:
        return all_valid, None


@app.command(name="version", help="Show version information")
def show_version() -> None:
    """Print version and exit."""
    console.print(
        Panel.fit(
            f"[bold cyan]py-gbcms[/bold cyan]\n"
            f"Version: [green]{__version__}[/green]\n"
            f"Python implementation of GetBaseCountsMultiSample (gbcms)",
            border_style="cyan",
            title="Version Info",
        )
    )
    raise typer.Exit()


def parse_bam_file(bam_string: str) -> tuple[str, str]:
    """
    Parse BAM file string in format SAMPLE:BAM_PATH.

    Args:
        bam_string: String in format "sample_name:bam_path"

    Returns:
        Tuple of (sample_name, bam_path)
    """
    parts = bam_string.split(":", 1)
    if len(parts) != 2:
        console.print(
            f"[red]Error:[/red] Incorrect format for --bam parameter: {bam_string}",
        )
        console.print("Expected format: SAMPLE_NAME:BAM_FILE")
        raise typer.Exit(1)
    return parts[0], parts[1]


def load_bam_fof(bam_fof_path: str) -> dict[str, str]:
    """
    Load BAM files from file-of-files.

    Args:
        bam_fof_path: Path to file containing sample names and BAM paths

    Returns:
        Dictionary mapping sample names to BAM paths
    """
    bam_files = {}
    with open(bam_fof_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                console.print(
                    f"[red]Error:[/red] Incorrect format at line {line_num} in {bam_fof_path}",
                )
                console.print("Expected format: SAMPLE_NAME<TAB>BAM_FILE")
                raise typer.Exit(1)

            sample_name, bam_path = parts
            if sample_name in bam_files:
                console.print(
                    f"[red]Error:[/red] Duplicate sample name: {sample_name}",
                )
                raise typer.Exit(1)

            bam_files[sample_name] = bam_path

    return bam_files


@count_app.command(name="run", help="Run base counting on variants")
def count_run(
    # Required arguments
    fasta: Annotated[
        Path,
        typer.Option(
            "--fasta",
            "-f",
            help="[bold cyan]Reference genome FASTA file[/bold cyan] (must be indexed with .fai)",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="📁 Required Input Files",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="[bold cyan]Output file path[/bold cyan]",
            rich_help_panel="📁 Required Input Files",
        ),
    ],
    # BAM input options
    bam: Annotated[
        list[str] | None,
        typer.Option(
            "--bam",
            "-b",
            help="BAM file in format [yellow]SAMPLE_NAME:BAM_FILE[/yellow] (can be specified multiple times)",
            rich_help_panel="🧬 BAM Input",
        ),
    ] = None,
    bam_fof: Annotated[
        Path | None,
        typer.Option(
            "--bam-fof",
            help="File containing sample names and BAM paths (tab-separated)",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="🧬 BAM Input",
        ),
    ] = None,
    # Variant input options (mutually exclusive)
    maf: Annotated[
        list[Path] | None,
        typer.Option(
            "--maf",
            help="Input variant file in [green]MAF format[/green] (can be specified multiple times)",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="🔬 Variant Input",
        ),
    ] = None,
    vcf: Annotated[
        list[Path] | None,
        typer.Option(
            "--vcf",
            help="Input variant file in [green]VCF format[/green] (can be specified multiple times)",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="🔬 Variant Input",
        ),
    ] = None,
    # Output format
    omaf: Annotated[
        bool,
        typer.Option(
            "--omaf",
            help="Output in MAF format (only with MAF input)",
            rich_help_panel="📤 Output Options",
        ),
    ] = False,
    positive_count: Annotated[
        bool,
        typer.Option(
            "--positive-count/--no-positive-count",
            help="Output positive strand counts (DPP/RDP/ADP)",
            rich_help_panel="📤 Output Options",
        ),
    ] = True,
    negative_count: Annotated[
        bool,
        typer.Option(
            "--negative-count/--no-negative-count",
            help="Output negative strand counts (DPN/RDN/ADN)",
            rich_help_panel="📤 Output Options",
        ),
    ] = False,
    fragment_count: Annotated[
        bool,
        typer.Option(
            "--fragment-count/--no-fragment-count",
            help="Output fragment counts (DPF/RDF/ADF)",
            rich_help_panel="📤 Output Options",
        ),
    ] = False,
    fragment_fractional_weight: Annotated[
        bool,
        typer.Option(
            "--fragment-fractional-weight",
            help="Use fractional weight (0.5) for fragments with disagreement",
            rich_help_panel="📤 Output Options",
        ),
    ] = False,
    # Quality filters
    maq: Annotated[
        int,
        typer.Option(
            "--maq",
            help="Mapping quality threshold",
            min=0,
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = 20,
    baq: Annotated[
        int,
        typer.Option(
            "--baq",
            help="Base quality threshold",
            min=0,
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = 0,
    filter_duplicate: Annotated[
        bool,
        typer.Option(
            "--filter-duplicate/--no-filter-duplicate",
            help="Filter reads marked as duplicate",
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = True,
    filter_improper_pair: Annotated[
        bool,
        typer.Option(
            "--filter-improper-pair/--no-filter-improper-pair",
            help="Filter reads marked as improperly paired",
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = False,
    filter_qc_failed: Annotated[
        bool,
        typer.Option(
            "--filter-qc-failed/--no-filter-qc-failed",
            help="Filter reads marked as QC failed",
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = False,
    filter_indel: Annotated[
        bool,
        typer.Option(
            "--filter-indel/--no-filter-indel",
            help="Filter reads containing indels",
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = False,
    filter_non_primary: Annotated[
        bool,
        typer.Option(
            "--filter-non-primary/--no-filter-non-primary",
            help="Filter non-primary alignments",
            rich_help_panel="🔍 Quality Filters",
        ),
    ] = False,
    # Performance options
    thread: Annotated[
        int,
        typer.Option(
            "--thread",
            "-t",
            help="Number of threads for parallel processing",
            min=1,
            rich_help_panel="⚡ Performance",
        ),
    ] = 1,
    backend: Annotated[
        str,
        typer.Option(
            "--backend",
            help="Parallelization backend: 'joblib' (default), 'loky', 'threading', or 'multiprocessing'",
            rich_help_panel="⚡ Performance",
        ),
    ] = "joblib",
    max_block_size: Annotated[
        int,
        typer.Option(
            "--max-block-size",
            help="Maximum number of variants per block",
            min=1,
            rich_help_panel="⚡ Performance",
        ),
    ] = 10000,
    max_block_dist: Annotated[
        int,
        typer.Option(
            "--max-block-dist",
            help="Maximum block distance in base pairs",
            min=1,
            rich_help_panel="⚡ Performance",
        ),
    ] = 100000,
    # Advanced options
    generic_counting: Annotated[
        bool,
        typer.Option(
            "--generic-counting",
            help="Use generic counting algorithm for complex variants",
            rich_help_panel="🔧 Advanced",
        ),
    ] = False,
    suppress_warning: Annotated[
        int,
        typer.Option(
            "--suppress-warning",
            help="Maximum number of warnings per type",
            min=0,
            rich_help_panel="🔧 Advanced",
        ),
    ] = 3,
    # Other options
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose",
            "-v",
            help="Enable verbose logging",
            rich_help_panel="🔧 Advanced",
        ),
    ] = False,
) -> None:
    """
    Calculate base counts in multiple BAM files for variants in VCF/MAF files.

    This tool counts the number of reference and alternate alleles at each variant
    position across multiple BAM files, with support for various quality filters
    and output formats.
    """
    # Setup logging
    # setup_logging(verbose)

    # Print banner
    console.print(
        Panel.fit(
            f"[bold cyan]py-gbcms[/bold cyan] v{__version__}\n"
            "Python implementation of GetBaseCountsMultiSample",
            border_style="cyan",
        )
    )

    # Validate inputs
    if not bam and not bam_fof:
        console.print(
            "[red]Error:[/red] Please specify at least one BAM file with --bam or --bam-fof",
        )
        raise typer.Exit(1)

    if not maf and not vcf:
        console.print(
            "[red]Error:[/red] Please specify at least one variant file with --maf or --vcf",
        )
        raise typer.Exit(1)

    if maf and vcf:
        console.print(
            "[red]Error:[/red] --maf and --vcf are mutually exclusive",
        )
        raise typer.Exit(1)

    # Parse BAM files
    bam_files = {}

    if bam:
        for bam_string in bam:
            sample_name, bam_path = parse_bam_file(bam_string)
            if sample_name in bam_files:
                console.print(
                    f"[red]Error:[/red] Duplicate sample name: {sample_name}",
                )
                raise typer.Exit(1)
            bam_files[sample_name] = bam_path

    if bam_fof:
        fof_bams = load_bam_fof(str(bam_fof))
        for sample_name, bam_path in fof_bams.items():
            if sample_name in bam_files:
                console.print(
                    f"[red]Error:[/red] Duplicate sample name: {sample_name}",
                )
                raise typer.Exit(1)
            bam_files[sample_name] = bam_path

    # Parse variant files
    variant_files = []
    input_is_maf = False
    input_is_vcf = False

    if maf:
        variant_files = [str(f) for f in maf]
        input_is_maf = True

    if vcf:
        variant_files = [str(f) for f in vcf]
        input_is_vcf = True

    # Validate input files before processing
    if not validate_input_files(fasta, bam_files, variant_files, input_is_maf, input_is_vcf)[0]:
        raise typer.Exit(1)

    # Display configuration
    config_table = Table(title="Configuration", show_header=False, border_style="cyan")
    config_table.add_column("Parameter", style="cyan")
    config_table.add_column("Value", style="green")

    config_table.add_row("Reference FASTA", str(fasta))
    config_table.add_row("Number of BAM files", str(len(bam_files)))
    config_table.add_row("Number of variant files", str(len(variant_files)))
    config_table.add_row("Input format", "MAF" if input_is_maf else "VCF")
    config_table.add_row("Output file", str(output))
    config_table.add_row("Threads", str(thread))
    config_table.add_row("Backend", backend)
    config_table.add_row("Mapping quality threshold", str(maq))
    config_table.add_row("Base quality threshold", str(baq))

    console.print(config_table)
    console.print()

    try:
        # Create configuration using legacy Config format (processor expects this)
        config = Config(
            fasta_file=str(fasta),
            bam_files=bam_files,
            variant_files=variant_files,
            output_file=str(output),
            mapping_quality_threshold=maq,
            base_quality_threshold=baq,
            filter_duplicate=filter_duplicate,
            filter_improper_pair=filter_improper_pair,
            filter_qc_failed=filter_qc_failed,
            filter_indel=filter_indel,
            filter_non_primary=filter_non_primary,
            output_positive_count=positive_count,
            output_negative_count=negative_count,
            output_fragment_count=fragment_count,
            fragment_fractional_weight=fragment_fractional_weight,
            max_block_size=max_block_size,
            max_block_dist=max_block_dist,
            num_threads=thread,
            backend=backend,
            input_is_maf=input_is_maf,
            input_is_vcf=input_is_vcf,
            output_maf=omaf,
            generic_counting=generic_counting,
            max_warning_per_type=suppress_warning,
        )

        # Process variants
        processor = VariantProcessor(config)
        processor.process()

        # Success message
        console.print()
        console.print(
            Panel.fit(
                "[bold green]✓[/bold green] Processing completed successfully!",
                border_style="green",
            )
        )

    except Exception as e:
        console.print()
        console.print(
            Panel.fit(
                f"[bold red]✗[/bold red] Error: {str(e)}",
                border_style="red",
            )
        )
        if verbose:
            console.print_exception()
        raise typer.Exit(1) from e


@validate_app.command(name="files", help="Validate input files")
def validate_files(
    fasta: Annotated[
        Path | None,
        typer.Option(
            "--fasta",
            "-f",
            help="Reference FASTA file to validate",
            rich_help_panel="Files to Validate",
        ),
    ] = None,
    bam: Annotated[
        list[str] | None,
        typer.Option(
            "--bam",
            "-b",
            help="BAM files to validate (SAMPLE:PATH format)",
            rich_help_panel="Files to Validate",
        ),
    ] = None,
    vcf: Annotated[
        list[Path] | None,
        typer.Option(
            "--vcf",
            help="VCF files to validate",
            rich_help_panel="Files to Validate",
        ),
    ] = None,
    maf: Annotated[
        list[Path] | None,
        typer.Option(
            "--maf",
            help="MAF files to validate",
            rich_help_panel="Files to Validate",
        ),
    ] = None,
) -> None:
    """
    Validate input files for gbcms.

    Checks:
    - File existence
    - Required indices (.fai for FASTA, .bai for BAM)
    - File format validity
    - Chromosome name consistency
    """

    console.print(
        Panel.fit(
            "[bold cyan]File Validation[/bold cyan]\n" "Checking input files for gbcms",
            border_style="cyan",
        )
    )

    results = Table(title="Validation Results", show_header=True, header_style="bold cyan")
    results.add_column("File Type", style="cyan")
    results.add_column("File Path", style="white")
    results.add_column("Status", style="white")
    results.add_column("Details", style="yellow")

    # Parse BAM files if provided
    bam_files = {}
    if bam:
        for bam_string in bam:
            sample_name, bam_path = parse_bam_file(bam_string)
            bam_files[sample_name] = bam_path

    # Parse variant files
    variant_files = []
    input_is_maf = False
    input_is_vcf = False

    if maf:
        variant_files = [str(f) for f in maf]
        input_is_maf = True

    if vcf:
        variant_files = [str(f) for f in vcf]
        input_is_vcf = True

    # Use the unified validation function with rich output
    # Note: fasta can be None, but if it is, the validation will handle it appropriately
    is_valid, results_table = validate_input_files(
        fasta or Path(""), bam_files, variant_files, input_is_maf, input_is_vcf, rich_output=True
    )

    # Handle results based on validation outcome
    if is_valid:
        console.print(
            Panel.fit(
                "[bold green]✓[/bold green] All files validated successfully!",
                border_style="green",
            )
        )
        raise typer.Exit(0)
    else:
        console.print(
            Panel.fit(
                "[bold red]✗[/bold red] Some files failed validation",
                border_style="red",
            )
        )
        raise typer.Exit(1)


@app.command(name="info", help="Show information about gbcms")
def show_info() -> None:
    """Display information about gbcms capabilities."""
    info_table = Table(title="gbcms Information", show_header=False, border_style="cyan")
    info_table.add_column("Category", style="bold cyan")
    info_table.add_column("Details", style="white")

    info_table.add_row("Version", __version__)
    info_table.add_row("Supported Input", "VCF, MAF")
    info_table.add_row("Supported Output", "VCF-like, MAF, Fillout")
    info_table.add_row("Variant Types", "SNP, DNP, Insertion, Deletion")
    info_table.add_row(
        "Quality Filters", "Mapping quality, Base quality, Duplicates, QC failed, etc."
    )
    info_table.add_row("Counting Methods", "DMP (default), Generic")
    info_table.add_row("Parallelization", "Multi-threaded with configurable threads")
    info_table.add_row("Dependencies", "pysam, numpy, typer, rich")

    console.print(info_table)
    console.print()

    console.print("[bold cyan]Example Usage:[/bold cyan]")
    console.print(
        "  gbcms count run --fasta ref.fa --bam s1:s1.bam --vcf vars.vcf --output out.txt"
    )
    console.print("  gbcms validate files --fasta ref.fa --bam s1:s1.bam")
    console.print("  gbcms version")


if __name__ == "__main__":
    app()
