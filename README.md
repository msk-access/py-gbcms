# py-gbcms - Python Implementation of gbcms

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-AGPL%203.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)

A high-performance Python reimplementation of [GetBaseCountsMultiSample](https://github.com/msk-access/GetBaseCountsMultiSample) for calculating base counts in multiple BAM files at variant positions specified in VCF or MAF files.

**Package**: `py-gbcms` | **CLI**: `gbcms`

## ✨ Features

- 🚀 **High Performance**: Multi-threaded processing with efficient BAM file handling
  - **Numba JIT compilation** for 50-100x speedup on counting operations
  - **joblib** for efficient local parallelization
  - **Ray** support for distributed computing across clusters
- 🎨 **Beautiful CLI**: Rich terminal output with progress bars and colored logging
- 🔒 **Type Safety**: Pydantic models for runtime validation and type checking
- 🐳 **Docker Support**: Containerized deployment for reproducibility
- 📊 **Multiple Formats**: Support for both VCF and MAF input/output formats
- 🧪 **Well Tested**: Comprehensive unit tests with high coverage
- 🔧 **Modern Python**: Built with type hints, Pydantic models, and modern Python practices

## 🚀 Quick Start

### Installation

```bash
# Install with all features
uv pip install "py-gbcms[all]"

# Or with pip
pip install "py-gbcms[all]"
```
**Requirements:** Python 3.11 or later

### Basic Usage

```bash
# Run
gbcms count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 8

# Run with MAF file
gbcms count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --maf variants.maf \
    --output counts.txt

### Docker Usage

```bash
docker pull ghcr.io/msk-access/getbasecounts:latest

# Run the container
docker run --rm \
    -v $(pwd)/data:/data \
    ghcr.io/msk-access/getbasecounts:latest \
    gbcms count run \
    --omaf \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.maf
```

### BAM File of Files Format

Create a tab-separated file (`bam_files.txt`):
```
sample1	/path/to/sample1.bam
sample2	/path/to/sample2.bam
sample3	/path/to/sample3.bam
```

Then use:

```bash
gbcms count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --vcf variants.vcf \

## Command Line Options

### Commands

gbcms uses subcommands for different operations:

- `gbcms count run`: Run base counting on variants (main command)
- `gbcms validate files`: Validate input files before processing
- `gbcms version`: Show version information
- `gbcms info`: Show tool capabilities and information

### Count Run Options

#### Required Arguments

- `--fasta`, `-f`: Reference genome FASTA file (must be indexed with .fai)
- `--output`, `-o`: Output file path

#### BAM Input (at least one required)

- `--bam`, `-b`: BAM file in format `SAMPLE_NAME:BAM_FILE` (can be specified multiple times)
- `--bam-fof`: File containing sample names and BAM paths (tab-separated)

#### Variant Input (one required)

- `--maf`: Input variant file in MAF format (can be specified multiple times)
- `--vcf`: Input variant file in VCF format (can be specified multiple times)

#### Output Options

- `--omaf`: Output in MAF format (only with MAF input)
- `--positive-count/--no-positive-count`: Output positive strand counts (default: True)
- `--negative-count/--no-negative-count`: Output negative strand counts (default: False)
- `--fragment-count/--no-fragment-count`: Output fragment counts (default: False)
- `--fragment-fractional-weight`: Use fractional weights for fragments (default: False)

#### Quality Filters
- `--maq`: Mapping quality threshold (default: 20)
- `--baq`: Base quality threshold (default: 0)
- `--filter-duplicate`: Filter duplicate reads (default: True)
- `--filter-improper-pair`: Filter improperly paired reads (default: False)
- `--filter-qc-failed`: Filter QC failed reads (default: False)
- `--filter-indel/--no-filter-indel`: Filter reads with indels (default: False)
- `--filter-non-primary/--no-filter-non-primary`: Filter non-primary alignments (default: False)

#### Performance Options
- `--thread`, `-t`: Number of threads (default: 1)
- `--max-block-size`: Maximum variants per block (default: 10000)
- `--max-block-dist`: Maximum block distance in bp (default: 100000)

#### Advanced Options
- `--generic-counting`: Use generic counting algorithm for complex variants
- `--suppress-warning`: Maximum warnings per type (default: 3)
- `--verbose`, `-v`: Enable verbose logging

### Validate Files Options

- `--fasta`, `-f`: Reference FASTA file to validate
- `--bam`, `-b`: BAM files to validate (SAMPLE:PATH format, multiple allowed)
- `--vcf`: VCF files to validate (multiple allowed)
- `--maf`: MAF files to validate (multiple allowed)

## Output Format

### VCF Output

The output is a tab-separated file with the following columns for each sample:

- `t_depth` (DP): Total depth
- `t_ref_count` (RD): Reference allele depth
- `t_alt_count` (AD): Alternate allele depth
- `t_ref_count_forward` (DPP): Positive strand depth (if enabled)
- `t_ref_count_forward` (RDP): Positive strand reference depth (if enabled)
- `t_alt_count_forward` (ADP): Positive strand alternate depth (if enabled)

### MAF Output

When using `--omaf`, the output maintains the MAF format with updated count columns.

## Development

### Setup Development Environment

```bash
# Clone the repository
git clone https://github.com/msk-access/py-gbcms.git
cd py-gbcms

# Install with development dependencies
uv pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=gbcms --cov-report=html

# Run specific test file
pytest tests/test_counter.py -v
```

### Code Quality

```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/

# Type checking
mypy src/
```

### Building Docker Image

```bash
# Build the image
docker build -t gbcms:latest .

# Run tests in container
docker run --rm gbcms:latest pytest
```

## Performance Comparison

Compared to the original C++ implementation:

| Feature | C++ Version | Python (Basic) | Python (Optimized) |
|---------|-------------|----------------|-------------------|
| Speed | ~1x | ~0.8-1.2x | ~2-5x** |
| Memory | Baseline | ~1.2x | ~1.5x |
| Multi-threading | OpenMP | concurrent.futures | joblib/Ray |
| Dependencies | bamtools, zlib | pysam, numpy | +numba, joblib, ray |
| Scalability | Single machine | Single machine | Multi-node clusters |

*Performance varies based on workload and Python version. Python 3.11+ shows significant improvements.

**With Numba JIT compilation and optimized parallelization. See [Fast VCF Parsing](docs/CYVCF2_SUPPORT.md) for benchmarks.

## Architecture

The package is organized into distinct modules with clear responsibilities:

```
py-gbcms/
├── src/gbcms/
│   ├── __init__.py
│   ├── cli.py              # 🎨 Typer CLI interface with Rich
│   ├── config.py           # ⚙️  Configuration dataclasses (legacy)
│   ├── models.py           # 🔒 Pydantic models for type safety ⭐
│   ├── variant.py          # 📄 Variant loading (VCF/MAF)
│   ├── counter.py          # 🐢 Pure Python counting (baseline)
│   ├── numba_counter.py    # ⚡ Numba-optimized counting (50-100x faster) ⭐
│   ├── parallel.py         # 🔄 joblib/Ray parallelization ⭐
│   ├── reference.py        # 🧬 Reference sequence handling
│   ├── output.py           # 📤 Output formatting
│   └── processor.py        # 🎯 Main processing pipeline
├── tests/                  # 🧪 Comprehensive test suite
├── scripts/                # 🛠️  Setup and verification scripts
├── Dockerfile              # 🐳 Production container
├── pyproject.toml          # 📦 Package configuration
└── docs/
    ├── README.md           # Main documentation
    ├── ARCHITECTURE.md     # Module relationships ⭐
    ├── INSTALLATION.md     # Setup guide
    ├── QUICKSTART.md       # 5-minute start
    ├── ADVANCED_FEATURES.md # Pydantic, Numba, Ray ⭐
    └── CLI_FEATURES.md     # CLI documentation
```

### Module Relationships

```
CLI (cli.py)
    ↓
Configuration (models.py)
    ↓
Processor (processor.py)
    ├─→ Variant Loader (variant.py)
    ├─→ Reference (reference.py)
    ├─→ Counting Engine
    │   ├─→ counter.py (Pure Python - flexible, slower)
    │   └─→ numba_counter.py (JIT compiled - 50-100x faster) ⭐
    ├─→ Parallelization (parallel.py)
    └─→ Output (output.py)
```

**Key Distinction:**
- **`counter.py`**: Pure Python, easy to debug, baseline performance
- **`numba_counter.py`**: JIT-compiled, 50-100x faster, for production

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for detailed module relationships.

## Documentation

📚 **[Complete Documentation](docs/README.md)** | [Quick Start](docs/QUICKSTART.md) | [Contributing Guide](CONTRIBUTING.md) | [Package Structure](docs/PACKAGE_STRUCTURE.md) | [Testing Guide](docs/TESTING_GUIDE.md)

- **User Guide**
  - [Input & Output Formats](docs/INPUT_OUTPUT.md)

- **Advanced**
  - [Advanced Features](docs/ADVANCED_FEATURES.md)
  - [Fast VCF Parsing (cyvcf2)](docs/CYVCF2_SUPPORT.md)

- **Reference**
  - [Architecture](docs/ARCHITECTURE.md)
  - [C++ Comparison](docs/CPP_FEATURE_COMPARISON.md)
  - [FAQ](docs/FAQ.md)

- **Docker & Deployment**
  - [Docker Guide](docs/DOCKER_GUIDE.md)

## Advanced Features

GetBaseCounts includes several advanced features for performance and scalability:

### 🔒 Type Safety with Pydantic

```python
from gbcms.models import GetBaseCountsConfig

# Runtime validation of all inputs
config = GetBaseCountsConfig(
    fasta_file=Path("reference.fa"),
    bam_files=[...],
    variant_files=[...],
)
```

### ⚡ Performance with Numba

```python
from gbcms.numba_counter import count_snp_batch

# JIT-compiled counting (50-100x faster)
counts = count_snp_batch(bases, qualities, positions, ...)
```

### 🔄 Parallelization with joblib

```bash
# Use joblib for efficient local parallelization
gbcms count run --thread 16 --backend joblib ...
```

### 🌐 Distributed Computing with Ray

```bash
# Install with Ray support
uv pip install "gbcms[ray]"

# Use Ray for distributed processing
gbcms count run --thread 32 --backend ray --use-ray ...
```

See [ADVANCED_FEATURES.md](ADVANCED_FEATURES.md) for detailed documentation and benchmarks.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Run the test suite and linters
5. Submit a pull request

## Citation

If you use this tool in your research, please cite:

```
GetBaseCountsMultiSample: A tool for calculating base counts in multiple BAM files
MSK-ACCESS Team
https://github.com/msk-access/GetBaseCountsMultiSample
```

## License

AGPL-3.0 License - See [LICENSE](LICENSE) for details.

## Support

- 🐛 Report bugs: [GitHub Issues](https://github.com/msk-access/getbasecounts/issues)
- 💬 Ask questions: [GitHub Discussions](https://github.com/msk-access/getbasecounts/discussions)
- 📧 Email: access@mskcc.org

## Acknowledgments

This is a Python reimplementation of the original C++ tool developed by the MSK-ACCESS team. Special thanks to the original authors and contributors.
