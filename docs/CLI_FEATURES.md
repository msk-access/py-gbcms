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
  gbcms count run --fasta ref.fa --bam s1:s1.bam --vcf vars.vcf --output out.txt
  gbcms validate files --fasta ref.fa --bam s1:s1.bam
  gbcms version
```

### Processing with Progress Bar

```bash
$ gbcms count run --fasta ref.fa --bam s1:s1.bam --vcf vars.vcf --output out.txt
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
