# CLI Features - gbcms

gbcms leverages **Typer** and **Rich** to provide a modern, user-friendly command-line interface with advanced features.

## Command Structure

The CLI uses a hierarchical subcommand structure for better organization:

```
gbcms
├── count
│   └── run          # Main counting command
├── validate
│   └── files        # Validate input files
├── version          # Show version info
└── info             # Show tool capabilities
```

## Key Typer Features Implemented

### 1. **Annotated Types with Rich Help Panels**

Options are organized into logical groups using `rich_help_panel`:

```python
fasta: Annotated[
    Path,
    typer.Option(
        "--fasta", "-f",
        help="[bold cyan]Reference genome FASTA file[/bold cyan]",
        rich_help_panel="📁 Required Input Files",
    ),
]
```

**Help panels include:**
- 📁 Required Input Files
- 🧬 BAM Input
- 🔬 Variant Input
- 📤 Output Options
- 🔍 Quality Filters
- ⚡ Performance
- 🔧 Advanced

### 2. **Multiple Values Support**

Multiple BAM files and variant files can be specified:

```bash
# Multiple --bam options
gbcms count run \
    --bam sample1:s1.bam \
    --bam sample2:s2.bam \
    --bam sample3:s3.bam \
    --vcf variants.vcf \
    --output out.txt

# Multiple variant files
gbcms count run \
    --vcf variants1.vcf \
    --vcf variants2.vcf \
    --vcf variants3.vcf \
    --output out.txt
```

Implementation:
```python
bam: Annotated[
    Optional[List[str]],
    typer.Option("--bam", "-b", help="..."),
] = None
```

### 3. **Subcommands for Different Operations**

#### Main Counting Command
```bash
gbcms count run --fasta ref.fa --bam s1:s1.bam --vcf vars.vcf --output out.txt
```

#### File Validation
```bash
gbcms validate files --fasta ref.fa --bam s1:s1.bam --vcf vars.vcf
```

#### Version Information
```bash
gbcms version
```

#### Tool Information
```bash
gbcms info
```

### 4. **Boolean Flags with Toggle Options**

Using Typer's flag syntax for clear enable/disable:

```bash
# Enable/disable filters
--filter-duplicate / --no-filter-duplicate
--positive-count / --no-positive-count
--fragment-count / --no-fragment-count
```

Implementation:
```python
filter_duplicate: Annotated[
    bool,
    typer.Option(
        "--filter-duplicate/--no-filter-duplicate",
        help="Filter reads marked as duplicate",
    ),
] = True
```

### 5. **Rich Markup in Help Text**

Help text uses Rich markup for better readability:

```python
help="BAM file in format [yellow]SAMPLE_NAME:BAM_FILE[/yellow]"
help="Input variant file in [green]MAF format[/green]"
help="[bold cyan]Reference genome FASTA file[/bold cyan]"
```

### 6. **Short and Long Options**

Common options have both short and long forms:

```bash
-f, --fasta          # Reference FASTA
-b, --bam            # BAM file
-o, --output         # Output file
-t, --thread         # Number of threads
-v, --verbose        # Verbose logging
```

### 7. **Input Validation**

Built-in validation using Typer's parameters:

```python
typer.Option(
    exists=True,          # File must exist
    file_okay=True,       # Must be a file
    dir_okay=False,       # Not a directory
    readable=True,        # Must be readable
    min=1,                # Minimum value
)
```

### 8. **No Args Shows Help**

The main app is configured to show help when no arguments are provided:

```python
app = typer.Typer(
    no_args_is_help=True,
    rich_markup_mode="rich",
)
```



## CLI Options Reference

This section documents all available command-line options for gbcms counting operations.

### Counting Control Options

#### Strand-Aware Counting
- `--strand-count/--no-strand-count` (default: True)
  - **Description**: Enable strand-specific counting including forward/reverse orientation
  - **Fields added**: DP_FORWARD, RD_FORWARD, AD_FORWARD, DP_REVERSE, RD_REVERSE, AD_REVERSE
  - **Use case**: Detect strand bias and orientation effects

#### Fragment-Aware Counting
- `--fragment-count/--no-fragment-count` (default: False)
  - **Description**: Enable fragment-based counting for paired-end data
  - **Fields added**: DPF, RDF, ADF, RDF_FORWARD, RDF_REVERSE, ADF_FORWARD, ADF_REVERSE
  - **Use case**: Analyze fragment orientation and improve accuracy for PE data

#### Fragment Weighting
- `--fragment-fractional-weight` (default: False)
  - **Description**: Use fractional weights (0.5) for fragments with orientation disagreement
  - **Use case**: Handle ambiguous fragment orientations more accurately

### Quality Filtering Options

#### Read-Level Filters
- `--filter-duplicate/--no-filter-duplicate` (default: False)
  - **Description**: Filter duplicate reads (marked with 0x400 flag)
- `--filter-improper-pair/--no-filter-improper-pair` (default: False)
  - **Description**: Filter improperly paired reads
- `--filter-qc-failed/--no-filter-qc-failed` (default: False)
  - **Description**: Filter reads that failed quality control
- `--filter-non-primary/--no-filter-non-primary` (default: False)
  - **Description**: Filter secondary and supplementary alignments

#### Quality Thresholds
- `--mapping-quality-threshold <int>` (default: 20)
  - **Description**: Minimum mapping quality for reads to be considered
- `--base-quality-threshold <int>` (default: 20)
  - **Description**: Minimum base quality for counting

#### Indel Filtering
- `--filter-indel/--no-filter-indel` (default: False)
  - **Description**: Filter reads containing insertions or deletions

### Input/Output Options

#### File Inputs
- `--fasta <file>`: Reference genome FASTA file (required)
- `--bam <sample_name>:<file>`: BAM file with sample name (can be specified multiple times)
- `--vcf <file>`: Input variant file in VCF format (can be specified multiple times)

#### Output Control
- `--output <file>`: Output file path (required)
- `--fillout`: Output in fillout format (extended MAF with all samples)

#### Input Format Detection
- `--input-is-maf`: Treat input as MAF format instead of VCF

### Performance Options

#### Parallelization
- `--threads <int>` (default: 1)
  - **Description**: Number of parallel threads for processing
- `--backend <str>` (default: "joblib")
  - **Options**: joblib, loky, threading, multiprocessing

#### Processing Control
- `--max-block-size <int>` (default: 10000)
  - **Description**: Maximum number of variants per processing block
- `--max-block-dist <int>` (default: 100000)
  - **Description**: Maximum distance between variants in a block

#### Optimization
- `--numba/--no-numba` (default: True)
  - **Description**: Use Numba JIT compilation for performance optimization

### Usage Examples

#### Basic Strand-Aware Counting
```bash
gbcms count \
  --fasta reference.fa \
  --bam tumor:tumor.bam \
  --bam normal:normal.bam \
  --vcf variants.vcf \
  --output results.vcf \
  --strand-count
```

#### Comprehensive Analysis with Filtering
```bash
gbcms count \
  --fasta reference.fa \
  --bam sample:sample.bam \
  --vcf variants.vcf \
  --output results.vcf \
  --strand-count \
  --fragment-count \
  --filter-duplicate \
  --mapping-quality-threshold 30 \
  --threads 4
```

#### Fragment-Aware Analysis
```bash
gbcms count \
  --fasta reference.fa \
  --bam sample:sample.bam \
  --vcf variants.vcf \
  --output results.vcf \
  --fragment-count \
  --fragment-fractional-weight \
  --filter-improper-pair
```

### Option Categories Summary

| Category | Options | Purpose |
|----------|---------|---------|
| **Counting** | `--strand-count`, `--fragment-count` | Control what types of counts to generate |
| **Filtering** | `--filter-*`, `--*-threshold` | Quality control and read filtering |
| **I/O** | `--fasta`, `--bam`, `--vcf`, `--output` | Specify input and output files |
| **Performance** | `--threads`, `--backend`, `--numba` | Control processing speed and resources |
| **Advanced** | `--fragment-fractional-weight`, `--max-block-size` | Fine-tune analysis parameters |

## Rich Integration Features

### 1. **Colored Output**

- **Cyan** for headers and important info
- **Green** for success messages
- **Red** for errors
- **Yellow** for warnings

### 2. **Panels and Boxes**

```python
console.print(
    Panel.fit(
        "[bold cyan]gbcms[/bold cyan]\n"
        "Version: [green]2.0.0[/green]",
        border_style="cyan",
    )
)
```

### 3. **Tables**

Configuration display:
```python
config_table = Table(title="Configuration", border_style="cyan")
config_table.add_column("Parameter", style="cyan")
config_table.add_column("Value", style="green")
```

Validation results:
```python
results = Table(title="Validation Results", header_style="bold cyan")
results.add_column("File Type", style="cyan")
results.add_column("Status", style="white")
```

### 4. **Progress Bars**

```python
with Progress(
    SpinnerColumn(),
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
) as progress:
    task = progress.add_task("Processing...", total=len(blocks))
```

### 5. **Rich Logging**

```python
logging.basicConfig(
    handlers=[RichHandler(
        console=console,
        rich_tracebacks=True,
        show_path=False
    )]
)
```

## Example Usage

### Basic Command with Organized Help

```bash
$ gbcms count run --help
```

Shows help organized into panels:
- 📁 Required Input Files
- 🧬 BAM Input
- 🔬 Variant Input
- 📤 Output Options
- 🔍 Quality Filters
- ⚡ Performance
- 🔧 Advanced

### File Validation

```bash
$ gbcms validate files \
    --fasta reference.fa \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --vcf variants.vcf
```

Output:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                          File Validation                                ┃
┃                  Checking input files for gbcms                 ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

                          Validation Results
┏━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━┓
┃ File Type┃ File Path               ┃ Status ┃ Details               ┃
┡━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━┩
│ FASTA    │ reference.fa            │ ✅ PASS│ File and index found  │
│ BAM      │ tumor:tumor.bam         │ ✅ PASS│ File and index found  │
│ BAM      │ normal:normal.bam       │ ✅ PASS│ File and index found  │
│ VCF      │ variants.vcf            │ ✅ PASS│ File found            │
└──────────┴─────────────────────────┴────────┴───────────────────────┘

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                  ✓ All files validated successfully!                    ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

### Version Information

```bash
$ gbcms version
```

Output:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                              Version Info                               ┃
┃                                                                          ┃
┃                          gbcms                                  ┃
┃                        Version: 2.0.0                                   ┃
┃          Python implementation of gbcms              ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

### Tool Information

```bash
$ gbcms info
```

Output:
```
                    gbcms Information
┏━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Version             ┃ 2.0.0                                            ┃
┃ Supported Input     ┃ VCF, MAF                                         ┃
┃ Supported Output    ┃ VCF-like, MAF, Fillout                           ┃
┃ Variant Types       ┃ SNP, DNP, Insertion, Deletion                    ┃
┃ Quality Filters     ┃ Mapping quality, Base quality, Duplicates, ...   ┃
┃ Counting Methods    ┃ DMP (default), Generic                           ┃
┃ Parallelization     ┃ Multi-threaded with configurable threads         ┃
┃ Dependencies        ┃ pysam, numpy, typer, rich                        ┃
┗━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

Example Usage:
  # VCF input (default: VCF output)
  gbcms count run --fasta ref.fa --bam sample1:tumor.bam --vcf vars.vcf --output out.txt
  # VCF input with fillout format
  gbcms count run --fasta ref.fa --bam sample1:tumor.bam --vcf vars.vcf --fillout --output out.txt
  # MAF input (default: sample-agnostic MAF output)
  gbcms count run --fasta ref.fa --bam sample1:tumor.bam --maf vars.maf --output out.txt
  # Multiple samples with explicit naming
  gbcms count run --fasta ref.fa --bam tumor.bam normal.bam --sample-name tumor,normal --maf vars.maf --output out.txt
  gbcms validate files --fasta ref.fa --bam sample1:sample1.bam
  gbcms version
```

### Processing with Progress Bar

```bash
$ gbcms count run --fasta ref.fa --bam sample1:tumor.bam --vcf vars.vcf --output out.txt
```

Output:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                          gbcms v2.0.0                           ┃
┃                Calculate base counts in multiple BAM files              ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

                            Configuration
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Reference FASTA              ┃ reference.fa                           ┃
┃ Number of BAM files          ┃ 1                                      ┃
┃ Number of variant files      ┃ 1                                      ┃
┃ Input format                 ┃ VCF                                    ┃
┃ Output file                  ┃ counts.txt                             ┃
┃ Threads                      ┃ 1                                      ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

[INFO] Loading reference sequence: reference.fa
[INFO] Loading variants file: variants.vcf
[INFO] 1000 variants loaded from file: variants.vcf
[INFO] Sorting variants
[INFO] Indexing variants
[INFO] Processing BAM file: sample1.bam

⠋ Processing sample1... ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:05

[INFO] Writing output to: counts.txt
[INFO] Successfully wrote 1000 variants to output file
[INFO] Finished processing

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                  ✓ Processing completed successfully!                   ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

## Benefits

1. **User-Friendly**: Clear, organized help with visual grouping
2. **Type-Safe**: Full type hints with validation
3. **Flexible**: Multiple input methods (individual files or file-of-files)
4. **Informative**: Rich feedback during processing
5. **Professional**: Beautiful terminal output
6. **Discoverable**: Subcommands make features easy to find
7. **Validated**: File validation before processing saves time

## Implementation Highlights

### Typer Features Used

- ✅ `Annotated` types for clean parameter definitions
- ✅ `rich_help_panel` for organized help
- ✅ Multiple values with `List[T]`
- ✅ Subcommands with `Typer()` instances
- ✅ Boolean flags with toggle syntax
- ✅ Path validation with `exists`, `file_okay`, etc.
- ✅ Rich markup in help text
- ✅ Short and long option names
- ✅ `no_args_is_help` for better UX

### Rich Features Used

- ✅ `Console` for colored output
- ✅ `Panel` for boxed messages
- ✅ `Table` for structured data
- ✅ `Progress` with spinners and bars
- ✅ `RichHandler` for beautiful logs
- ✅ Rich markup in strings
- ✅ Exception formatting

This creates a modern, professional CLI that's both powerful and pleasant to use! 🎨
