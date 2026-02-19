# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.6.1] - 2026-02-19

### 🔧 Fixed
- **Per-haplotype trimming**: Fixed `slice index starts at 7 but ends at 6` panic in `counting.rs` on asymmetric indels. Replaced shared symmetric trim with independent per-haplotype `trim_haplotype()` function that calculates bounds safely for each allele

### ✨ Added
- **Tolerant REF validation**: Variants with ≥90% REF match against the FASTA are now counted (status `PASS_WARN_REF_CORRECTED`) instead of being silently rejected. The FASTA REF is used for haplotype construction. Variants with <90% match are still rejected as `REF_MISMATCH`

### 📚 Documentation
- **Visual posters**: Added overview, normalization, and read-filter/counting-metrics posters (JPG) to reference documentation pages with lightbox support
- **Embedded PDFs**: Added inline PDF viewer for allele classification guide and detailed overview presentation via `mkdocs-pdf` plugin
- **Variant normalization**: Updated REF validation docs with 3-tier flowchart, `PASS_WARN_REF_CORRECTED` status, and EGFR exon 19 real-world example

### 🔧 CI
- **`deploy-docs.yml`**: Added `mkdocs-pdf` to docs CI pip install dependencies

## [2.6.0] - 2026-02-18

### ✨ Added
- **Adaptive context padding**: Dynamically increases `ref_context` flanking in tandem repeat regions (homopolymer through hexanucleotide). Formula: `max(default, repeat_span/2 + 3)`, capped at 50bp. Enabled by default (`--adaptive-context/--no-adaptive-context`)
- **`gbcms normalize` command**: Standalone variant normalization (left-align + REF validate) without counting, outputs TSV with original and normalized coordinates
- **Nextflow parameters**: `fragment_qual_threshold`, `context_padding`, `show_normalization`, `adaptive_context` now configurable in `nextflow.config`
- **Docs restructure**: Split monolithic `variant-counting.md` into 4 focused pages: Variant Normalization, Allele Classification, Counting Metrics, Read Filters
- **HPC install docs**: Micromamba-based source install with Python 3.13

### 🔧 Fixed
- **Interior REF guard** for large deletions (>50bp): Reads falling entirely within a deleted region are now correctly classified as REF instead of ALT by Smith-Waterman
- **Windowed reciprocal overlap**: Improved shifted indel detection using bidirectional overlap scoring
- **Complex variant counting** (EPHA7 `TCC→CT`): Fixed base quality extraction for all variant-type handlers (`check_insertion`, `check_deletion`, `check_mnp`, `check_complex`)
- **MAF VCF-style conversion**: Corrected complex variant handling in MAF→internal coordinate conversion
- **Lint**: Fixed ruff I001/E402/B905 and black formatting in `pipeline.py`

### 🔄 Changed
- **Dead code removed**: `GenomicInterval` class, `Variant.interval` property, `fragment_counting` config field
- **`fetch_single_base()` refactored**: Delegates to `fetch_region()`, removing 33 lines of duplicated chr-prefix retry logic
- **Release guide**: Updated version locations table with exact line numbers and verification command
- **Nextflow pipeline diagram** added to docs index

## [2.5.0] - 2026-02-12

### ✨ Added
- **`--preserve-barcode` flag**: Keeps original `Tumor_Sample_Barcode` from input MAF instead of overriding with BAM sample name (MAF→MAF workflows)
- **`--column-prefix` parameter**: Controls prefix for gbcms count columns in MAF output (default: none; use `--column-prefix t_` for legacy compatibility) ⚠️
- **`CoordinateKernel`**: Centralized MAF↔internal 0-based coordinate conversion with variant-type-aware logic for SNP, insertion, deletion, and complex variants
- **Nextflow `FILTER_MAF` module**: Per-sample MAF variant filtering by `Tumor_Sample_Barcode` supporting exact match, regex, and multi-select (comma-separated) modes
- **Nextflow `PIPELINE_SUMMARY` module**: Aggregated per-sample filtering statistics with formatted console output
- **Nextflow `--filter_by_sample` parameter** and samplesheet `tsb` column for multi-sample MAF workflows
- **Nextflow documentation**: Samplesheet `tsb` column guide, `--filter_by_sample` parameter reference, multi-sample MAF filtering examples

### 🔧 Fixed
- **Fragment quality extraction** (critical): All variant-type handlers (`check_insertion`, `check_deletion`, `check_mnp`, `check_complex`) now return actual `base_qual` from CIGAR walk instead of 0 — fixes systematic ALT undercount for indels in fragment-level consensus
- **FILTER_MAF heredoc conflict**: Restructured script from Python shebang to `python3 << 'PYEOF'` pattern, resolving `SyntaxError` from bash syntax in Python context
- **FILTER_MAF string quoting**: Changed to single-quoted Python strings for Nextflow variable interpolation to prevent CSV-parsed double-quote conflicts
- **`splitCsv` quote handling**: Added `quote:'"'` parameter for correct RFC 4180 parsing of comma-separated TSB values within quoted CSV fields
- **mypy `no-redef` error**: Removed redundant type annotation in `output.py` else branch

### 🔄 Changed
- **MAF output column prefix default**: Changed from `t_` to empty string (no prefix). Use `--column-prefix t_` for legacy `t_ref_count` / `t_alt_count` style columns ⚠️
- **`MafWriter` refactored**: MAF→MAF path preserves all original columns verbatim; VCF→MAF path builds row from GDC fieldnames with `CoordinateKernel` coordinate conversion
- **Nextflow `GBCMS_RUN` input**: Variants bundled into sample tuple `(meta, bam, bai, variants)` instead of separate channel
- **Nextflow `GBCMS` workflow**: Simplified to 2-channel interface (`ch_samples`, `ch_fasta`) from 3 channels
- **Nextflow config**: Added `column_prefix`, `preserve_barcode`, `filter_by_sample` parameters

## [2.4.0] - 2026-02-10

### ✨ Added
- **Fragment Consensus Engine**: `FragmentEvidence` struct with u64 QNAME hashing and quality-weighted R1/R2 consensus; ambiguous fragments are discarded (not assigned to REF)
- **`--fragment-qual-threshold`**: New CLI option (default 10) controlling consensus quality difference for fragment conflict resolution
- **Windowed Indel Detection**: ±5bp positional scan with 3-layer safeguards (sequence identity, closest match, reference context validation)
- **Quality-Aware Complex Matching**: Masked comparison that ignores bases below `--min-baseq`; 3-case comparison (equal-length, ALT-only, REF-only) with ambiguity detection
- **Variant Counting Guide**: New `docs/reference/variant-counting.md` with algorithm diagrams for all variant types (~700 lines)
- **MAF Normalization Docs**: Added indel normalization and coordinate handling to `docs/reference/input-formats.md`
- **47 Tests**: Up from 16 — added `test_shifted_indels.py` (15), `test_fuzzy_complex.py` (14), `test_fragment_consensus.py` (4)

### 🔄 Changed
- **`--min-baseq` default**: `0` → `20` (Phred Q20) — activates quality masking by default for improved accuracy on low-coverage samples ⚠️
- **`--version` flag**: Added to CLI (`gbcms --version`)
- **Deploy-Docs Workflow**: Replaced `mkdocs gh-deploy` with `mike` for multi-version documentation; deploys `stable` (tagged version) from main and `dev` from develop branch; added `extra.version.provider: mike` to `mkdocs.yml` with version switcher widget

### 🔧 Fixed
- **Fragment double-counting bug**: R1+R2 pairs previously counted as two independent observations; now collapsed via quality-weighted consensus
- **MAF Input Hardening**: Graceful handling of missing/malformed fields with warnings instead of crashes
- **CI Release Pipeline**: Stabilized manylinux builds — migrated from manylinux_2_28 to manylinux_2_34, resolved OpenSSL/CURL vendor conflicts via `docker-options` pattern
- **Type stubs**: `_rs.pyi` and `gbcms_rs.pyi` synced with Rust bindings (added `ref_context`, `ref_context_start`, `fragment_qual_threshold`)
- **Linting**: All files pass black, ruff, and mypy

### 📚 Documentation
- Architecture comparison table updated with windowed indels and masked comparison
- Nextflow config and docs updated with new `min_baseq` default
- Testing guide expanded with Phase 2a/2b test files
- HPC/RHEL 8 installation instructions updated with `clangdev` header management
- Release guide updated with docs version locations
- `.antigravity` project files updated with current Rust LOC (~1270) and test counts (47)

## [2.3.0] - 2026-02-06

### ✨ Added
- **Nextflow BAI Auto-Discovery**: Checks `.bam.bai` and `.bai` extensions automatically
- **Documentation Modernization**: Hierarchical navigation, glightbox, panzoom, abbreviations
- **Performance Benchmarks**: cfDNA duplex sample metrics in documentation
- **RHEL 8 Installation Guide**: Conda-based source installation for legacy Linux

### 🔄 Changed
- **Dockerfile**: Added `procps`, `bash`, OCI labels, `maturin[patchelf]`, selective COPY
- **Nextflow Config**: `--platform linux/amd64`, shell config, local profile, observability (trace/report/timeline/dag)
- **MkDocs**: Switched to `navigation.sections`, 20+ abbreviations with hover tooltips
- **GitHub Actions**: Consolidated deploy-docs workflows, added caching and PR validation
- **CI Wheels**: Migrated from `manylinux_2_28` to `manylinux_2_34` (AlmaLinux 9 with OpenSSL 3.0+)

### 🔧 Fixed
- **Nextflow**: Empty `--suffix` argument no longer causes failures
- **Admonitions**: Converted GitHub-style alerts to MkDocs syntax
- **CI Build**: Resolved `curl-sys` OpenSSL version conflict by switching to manylinux_2_34

## [2.2.0] - 2026-02-04

### ✨ Added
- **Multi-platform Wheel Publishing**: Maturin-based CI builds for Linux (x86_64, aarch64), macOS (Intel, Apple Silicon), and Windows
- **Structured Logging**: New `utils/logging.py` module with Rich console output, timing utilities, and log file support
- **Mermaid Diagrams**: Architecture documentation with interactive flowcharts
- **Release Guide**: Comprehensive `docs/RELEASE.md` with git-flow workflow

### 🔄 Changed
- **Folder Restructure**: Moved Rust code to `rust/` (bundled as `gbcms._rs`)
- **Config Hierarchy**: Nested Pydantic models (`ReadFilters`, `QualityThresholds`, `OutputConfig`) for better organization
- **Code Quality**: Added `__all__` exports, docstrings, and type hints across all modules
- **StrEnum**: Modern enum pattern with Python 3.10 backport

### 📚 Documentation
- New `docs/ARCHITECTURE.md` with system diagrams
- New `docs/DEVELOPMENT.md` (developer guide)
- New `docs/TESTING.md` (testing guide)
- Updated MkDocs with mermaid2 plugin and snippet includes

## [2.1.2] - 2025-11-25

### 🔧 Fixed
- **PyPI Distribution**: Fixed source distribution size issue by correctly excluding large files (tests, docs, etc.) via `pyproject.toml` configuration.

## [2.1.1] - 2025-11-25 [YANKED]

!!! warning "Yanked Release"
    This release was yanked from PyPI due to a source distribution size limit error. Use 2.1.2 instead.

### 🔧 Fixed
- **PyPI Distribution**: Added MANIFEST.in (failed to work with Hatchling) to reduce source distribution size
- **Documentation**: Added comprehensive Installation guide
- **Documentation**: Unified Contributing guide (merged code + docs contributions)
- **Documentation**: Added Changelog to documentation navigation

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
- `gbcms._rs`: PyO3-based Rust extension (bundled in wheel)
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
