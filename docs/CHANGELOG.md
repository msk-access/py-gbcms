# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2025-11-25

### ✨ Added

#### Nextflow Workflow
- **Production-ready Nextflow workflow** for processing multiple samples in parallel
- **SLURM cluster support** with customizable queue configuration
- **Per-sample suffix support** via optional `suffix` column in samplesheet
- **Docker and Singularity profiles** for containerized execution
- **Automatic BAI index discovery** with validation
- **Resume capability** for failed workflow runs
- **Resource management** with automatic retry and scaling
- **Comprehensive documentation** in `docs/NEXTFLOW.md` and `nextflow/README.md`

#### Documentation
- **Usage pattern comparison** guide (`docs/WORKFLOWS.md`) for choosing between CLI and Nextflow
- **MkDocs integration** for beautiful GitHub Pages documentation
- **Local documentation preview** with live reload (`mkdocs serve`)
- **Staging deployment** from `develop` branch for testing docs
- **Production deployment** from `main` branch
- **Reorganized documentation structure** with clear CLI vs Nextflow separation
- **CLI Quick Start guide** (`docs/quick-start.md`)

### 🔧 Changed
- **Documentation workflow**: docs now live on `main` branch with automated deployment
- **GitBook integration**: configured to read from `main` branch
- **Nextflow module**: improved parameter passing with meta.suffix support

### 📝 Documentation
- Complete Nextflow workflow guide with SLURM examples
- Per-sample suffix usage examples
- Git-flow documentation workflow guide
- Local preview instructions
- Updated README with clear usage pattern separation

## [2.0.0] - 2025-11-21

### 🚀 Major Rewrite

Version 2.0.0 represents a complete rewrite of py-gbcms with a focus on performance, correctness, and modern architecture.

### ✨ Added

#### Core Features
- **Rust-based Counting Engine**: Hybrid Python/Rust architecture for 20x+ performance improvement
- **Strand Bias Statistics**: Fisher's exact test p-values and odds ratios for both reads (`SB_PVAL`, `SB_OR`) and fragments (`FSB_PVAL`, `FSB_OR`)
- **Fragment-Level Counting**: Majority-rule fragment counting with strand-specific counts (`RDF`, `ADF`)
- **Variant Allele Fractions**: Read-level (`VAF`) and fragment-level (`FAF`) allele fraction calculations
- **Thread Control**: Explicit control over parallelism via `--threads` argument (default: 1)

#### Input/Output
- **VCF Output Format**: Standard VCF with comprehensive INFO and FORMAT fields
- **MAF Output Format**: Extended MAF with custom columns for strand counts and statistics
- **Column Preservation**: Input MAF columns are preserved in output
- **Multiple BAM Support**: Process multiple samples via `--bam-list` or repeated `--bam` arguments
- **Sample ID Override**: Explicit sample naming via `--bam sample_id:path` syntax

#### Filters
- `--filter-duplicates`: Filter duplicate reads (default: enabled)
- `--filter-secondary`: Filter secondary alignments
- `--filter-supplementary`: Filter supplementary alignments
- `--filter-qc-failed`: Filter reads that failed QC
- `--filter-improper-pair`: Filter improperly paired reads
- `--filter-indel`: Filter reads with indels in CIGAR

#### CLI & Usability
- **Modern CLI**: Built with Typer and Rich for beautiful terminal output
- **Progress Tracking**: Real-time progress bars and status indicators
- **Direct Invocation**: Use `gbcms run` instead of `python -m gbcms.cli`
- **Output Customization**: `--suffix` flag for output filename customization
- **Flexible Input**: Support for both VCF and MAF input formats

#### Infrastructure
- **Docker Support**: Production-ready multi-stage Dockerfile with optimized layers
- **Type Safety**: Full type annotations with mypy support
- **Type Stubs**: Provided `.pyi` stub file for Rust extension
- **Comprehensive Tests**: Extended test suite with accuracy and filter validation
- **CI/CD**: GitHub Actions workflows for testing, linting, and releases

### 🔄 Changed

#### Architecture
- Migrated from pure Python to hybrid Python/Rust architecture
- Core counting logic implemented in Rust using `rust-htslib`
- Data parallelism over variants with per-thread BAM readers

#### Output Formats
- **VCF FORMAT fields**: Strand-specific counts now use comma-separated values (e.g., `RD=5,3` for forward,reverse)
- **MAF columns**: Standardized column names (`t_ref_count_forward`, `t_alt_count_reverse`, etc.)
- **Coordinate System**: Internal 0-based indexing with correct conversion for VCF (1-based) and MAF output

#### Performance
- **Speed**: 20x+ faster than v1.x on typical datasets
- **Memory**: Efficient per-thread BAM readers with minimal overhead
- **Scalability**: Configurable thread pool for optimal resource usage

#### Dependencies
- **Python**: Updated to require Python ≥3.10
- **Rust**: pyo3 0.27.1, rust-htslib 0.51.0, statrs 0.18.0
- **Python Packages**: pysam ≥0.21.0, typer ≥0.9.0, rich ≥13.0.0, pydantic ≥2.0.0

### 🗑️ Removed

- **Legacy Python Counting**: Pure Python implementation removed in favor of Rust
- **Old CLI**: Deprecated `python -m gbcms.cli` entry point
- **Unused Dependencies**: Removed `cyvcf2` and `numba` (no longer needed)
- **Pre-commit Hooks**: Removed in favor of explicit linting in CI

### 🐛 Fixed

- Correct handling of complex variants (MNPs, DelIns)
- Proper strand assignment for fragment counting
- Reference validation against FASTA for all variant types
- Thread-safe BAM access with per-thread readers

### 📚 Documentation

- Complete rewrite of all documentation
- New guides: `INSTALLATION.md`, `CLI_FEATURES.md`, `INPUT_OUTPUT.md`
- Comprehensive API documentation
- Docker usage examples
- Contributing guidelines updated

### 🔧 Technical Details

#### Rust Components
- `gbcms_rs`: PyO3-based extension module
- Fisher's exact test via `statrs` crate
- Rayon-based parallelism with configurable thread pools
- Safe memory management with Rust's ownership model

#### Testing
- 16 comprehensive test cases
- Accuracy validation with synthetic BAM files
- Filter validation for all read flag combinations
- Integration tests with real-world data

### ⚠️ Breaking Changes

Version 2.0.0 is **not backward compatible** with 1.x. Key breaking changes:

1. **CLI syntax**: Use `gbcms run` instead of `python -m gbcms.cli`
2. **Output format**: VCF/MAF column structures have changed
3. **Default behavior**: Only duplicate filtering enabled by default (was: all filters)
4. **Dependencies**: Requires Rust toolchain for installation from source
5. **Python version**: Minimum Python 3.10 (was: 3.8)

### 📦 Installation

```bash
# From PyPI (includes pre-built wheels)
pip install py-gbcms

# From source (requires Rust)
pip install git+https://github.com/msk-access/py-gbcms.git

# Docker
docker pull ghcr.io/msk-access/py-gbcms:2.0.0
```

### 🙏 Acknowledgments

This rewrite was designed and implemented with a focus on correctness, performance, and modern best practices in bioinformatics software development.

---

## [1.x] - Legacy

Previous versions (1.x) used a pure Python implementation. See git history for details.
