# gbcms Development Guide

## Quick Start

```bash
# Clone + build (builds Rust extension in release mode)
git clone https://github.com/msk-access/gbcms.git && cd gbcms
python -m venv .venv && source .venv/bin/activate
maturin develop --release
gbcms --version
```

## Project Structure

```
gbcms/
├── src/gbcms/          # Python package
│   ├── cli.py          # CLI entry point — 4-layer validation
│   ├── pipeline.py     # Orchestration, Parquet dispatch
│   ├── normalize.py    # Standalone normalization workflow
│   ├── io/             # VcfReader, MafReader, VcfWriter, MafWriter
│   ├── core/           # CoordinateKernel (VCF/MAF → 0-based)
│   ├── models/         # GbcmsConfig, OutputConfig, AlignmentConfig
│   └── utils/          # Logging utilities
├── rust/               # Rust extension (gbcms._rs)
│   └── src/
│       ├── counting/   # engine, variant_checks, alignment, pairhmm,
│       │               #   fragment, mfsd, utils, mod
│       ├── normalize/  # engine, left_align, decomp, fasta, repeat,
│       │               #   types, mod
│       ├── parquet_writer.rs  # write_fsd_parquet() — Arrow+ZSTD
│       ├── types.rs    # Variant, BaseCounts (PyO3 bindings)
│       ├── stats.rs    # Fisher's exact test
│       └── lib.rs      # PyO3 module exports
├── tests/              # 15 test files — see ARCHITECTURE.md
├── docs/               # MkDocs documentation site
├── nextflow/           # Nextflow pipeline + modules
└── pyproject.toml      # maturin build config
```

## Build

```bash
maturin develop           # dev build (unoptimized, fast)
maturin develop --release # release build (always use before regression tests)
maturin build --release --out dist  # wheel for distribution
```

## Full Lint Suite (Run Before Every Commit)

```bash
source .venv/bin/activate
ruff check src/ tests/
black --check src/ tests/
mypy src/gbcms/ tests/
cargo clippy --manifest-path rust/Cargo.toml -- -D warnings
```

All four must be clean. Fix in this order: ruff → black → mypy → clippy.

**After mypy fixes in test files** using `pysam.AlignedSegment`:
- Use `.cigartuples = [(...)]` (list of tuples), NOT `.cigar = (...)`
- Add `# type: ignore[assignment]` to `.query_qualities = [...]` lines
  (pysam stubs say `array[Any]` but lists work at runtime)

## Testing

```bash
pytest --no-cov -q          # fast run (no coverage)
pytest -v                   # verbose
pytest tests/test_mfsd_flag.py -v  # specific file
cargo test --manifest-path rust/Cargo.toml  # Rust inline tests
```

### Key Test Invariants

Every counting test must assert:
```python
assert counts.dp >= counts.rd + counts.ad      # DP includes 'neither' reads
assert counts.dpf >= counts.rdf + counts.adf   # DPF includes discarded
assert counts.rd == counts.rd_fwd + counts.rd_rev  # strand consistency
assert counts.ad == counts.ad_fwd + counts.ad_rev
```

### mFSD Tests
```python
# Column counts: exact numbers matter
assert len([c for c in maf.columns if c.startswith("mfsd_")]) == 34
assert len(maf.columns) == 145  # no --mfsd
assert len(maf.columns) == 179  # with --mfsd
```

## Git Workflow (git-flow)

```
main          ← stable releases (tagged)
  └── develop      ← integration branch (target for PRs)
       └── feature/*   ← individual features/fixes
```

- Branch from `develop`: `git checkout -b feature/my-thing develop`
- Merge into `develop` with `--no-ff`; never push directly to `main`
- Releases: `release/X.Y.Z` branch → bump version → PR to `main` → tag

## CLI Validation Rules (4 Layers — Enforce This Order)

1. **Parse-time (Typer)**: `Enum` for constrained choices, `min=`/`max=` for ranges. Typer rejects before Python code runs.
2. **Pre-model (cli.py command body)**: File-extension checks, cross-option dependencies (e.g. `--mfsd-parquet` requires `--mfsd`). Log at `ERROR`, raise `typer.Exit(1)`.
3. **Model-time (Pydantic)**: `Field(ge=..., le=...)` and `@model_validator` in `models/core.py`.
4. **No silent skips**: Missing inputs → fail-fast or explicit opt-out (`--lenient-bam`). A missing required input silently continuing is a bug.

The `cli.py` module docstring documents this order. Keep it in sync.

## mFSD Implementation Rules

- `--mfsd` gates output columns/fields — writers check `self.mfsd`
- `--mfsd-parquet` always requires `--mfsd` — validated at CLI and Pydantic
- `ref_sizes` and `alt_sizes` in `BaseCounts` are **Rust-internal** (no `#[pyo3(get)]`)
- `write_fsd_parquet()` is called from `pipeline.py` after `count_bam()`, not from within the counting engine
- Parquet uses ZSTD(1) compression (`features = ["arrow", "zstd"]` in Cargo.toml) — Snappy/other codecs are disabled by `default-features = false`

## Code Standards

### Python
- All public functions: type hints + Google-style docstrings
- `__all__` in every module
- Use `logging` — never `print()` or `console.print()` for status
- Pydantic models for all config (no raw dicts passed between layers)
- `SimpleNamespace` only for test fixtures (`_zero_counts()`), never production

### Rust
- `///` doc comments on all public items
- `anyhow::Result` for propagated errors; `pyo3::exceptions::PyIOError` for errors crossing the FFI boundary
- `log` crate for logging (`debug!`, `info!`, `warn!`)
- No `unwrap()` in production code paths — use `?` or explicit error handling

## QC Checklist (Per Phase — Before Implementing, Review and Check)

### Logging
- [ ] All user-facing messages use proper log levels (`debug`/`info`/`warning`/`error`)
- [ ] No `print()` or `console.print()` for status
- [ ] Critical operations have timing/performance logging
- [ ] Errors include context (file path, sample name, position)

### Commenting
- [ ] All modules have module-level docstrings stating purpose and key design decisions
- [ ] All public functions have docstrings with Args/Returns
- [ ] Complex logic has inline comments explaining *why* (not what)
- [ ] No commented-out code

### No Code Duplication
- [ ] Common patterns extracted to helper functions in `utils/`
- [ ] Similar code across files consolidated
- [ ] No copy-pasted logic between MafWriter and VcfWriter

### No Unused Code
- [ ] All imports are used (ruff catches this)
- [ ] No dead functions/methods
- [ ] `__all__` exports match public API

### Tests
- [ ] New features have corresponding test file or added to existing test
- [ ] All 4 invariants asserted in counting tests
- [ ] `pytest --no-cov -q` passes (107+ tests)
- [ ] mFSD column counts verified if output.py changed

## Phase Implementation Workflow

```
1. Review all affected files against QC checklist above
2. Document findings in implementation plan (in brain/ artifacts)
3. Fix QC issues FIRST before adding new features
4. Implement phase changes
5. Run full lint suite (ruff, black, mypy, clippy)
6. Run tests: pytest --no-cov -q
7. Run integration test with 22-BAM real data (for counting engine changes)
8. Re-check QC standards
9. Commit with detailed message describing what/why
```

## Real-Data Integration Test (22-BAM Suite)

For any changes to counting, mFSD, or output writers:

```bash
# Build BAM list
for bam in ~/Downloads/bams_slice/P-*.bam; do
    name="$(basename "$bam" .bam)"; name="${name%_subset}"
    echo "$name $bam"
done > /tmp/gbcms_bam_list.txt

# Run without mFSD (regression check)
gbcms run \
  --variants ~/Downloads/bams_slice/filtered_indels.maf \
  --bam-list /tmp/gbcms_bam_list.txt \
  --fasta /path/to/Homo_sapiens_assembly19.fasta \
  --output-dir /tmp/gbcms_no_mfsd \
  --format maf --preserve-barcode --lenient-bam --threads 4

# Run with mFSD + Parquet
gbcms run ... --mfsd --mfsd-parquet --output-dir /tmp/gbcms_with_mfsd

# Verify: 22 samples, 0 failed, exit 0
# No --mfsd: 145 columns, 0 parquet files
# --mfsd: 179 columns (+34 mFSD), 22 .fsd.parquet files (62–208 KB each)
# Baseline regression: alt_count/ref_count bit-for-bit identical to previous run
```

## Environment Variables

```bash
GBCMS_LOG_LEVEL=DEBUG RUST_LOG=debug gbcms run ...
```

## Release Process

1. `git checkout -b release/X.Y.Z develop`
2. Bump `version` in `pyproject.toml` and `rust/Cargo.toml`
3. Update `CHANGELOG.md`
4. Commit, PR to `main`, merge, tag: `git tag X.Y.Z && git push origin X.Y.Z`
5. CI: build → PyPI publish → Docker push to GHCR
6. Merge tag back to `develop`
