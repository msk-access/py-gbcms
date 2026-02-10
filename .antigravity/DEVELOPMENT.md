# py-gbcms Development Guide

## Quick Start

```bash
# Install in editable mode (builds Rust extension)
maturin develop --release

# Run tests
pytest tests/

# Run CLI
gbcms run -v variants.vcf -b sample.bam -f reference.fa -o output/
```

## Project Structure

```
py-gbcms/
├── src/gbcms/          # Python package
│   ├── cli.py          # Entry point
│   ├── pipeline.py     # Orchestration
│   ├── io/             # I/O adapters
│   ├── core/           # Coordinate kernel
│   ├── models/         # Pydantic models
│   └── utils/          # Logging utilities
├── rust/               # Rust extension (top-level)
│   ├── src/
│   │   ├── counting.rs # BAM processing (~1070 LOC)
│   │   ├── types.rs    # PyO3 bindings
│   │   ├── stats.rs    # Fisher's test
│   │   └── lib.rs      # Module exports
│   └── Cargo.toml
├── tests/              # Test suite (47 tests)
│   ├── test_accuracy.py
│   ├── test_shifted_indels.py   # Windowed indel detection
│   ├── test_fuzzy_complex.py    # Masked complex matching
│   ├── test_filters.py
│   ├── test_strand_counts.py
│   └── ...
└── docs/               # Documentation
```

## Building (Unified Python + Rust)

```bash
# Development build (editable install with Rust)
maturin develop --release

# Build wheel
maturin build --release --out dist
```

## Code Quality Standards

- All modules have `__all__` exports
- All public functions have docstrings
- Use `logging` module (not print)
- Pydantic models for configuration
- Type hints throughout

---

## QC Review Checklist (Per Phase)

**Before implementing any phase, review the codebase to ensure:**

### 1. Logging
- [ ] All user-facing messages use proper logging levels (`debug`/`info`/`warning`/`error`)
- [ ] No `print()` statements or `console.print()` for status (use `logger`)
- [ ] Critical operations have timing/performance logging
- [ ] Errors include context for debugging

### 2. Commenting
- [ ] All modules have docstrings describing purpose
- [ ] All public functions have docstrings with Args/Returns
- [ ] Complex logic has inline comments explaining "why"
- [ ] No commented-out code (delete or use version control)

### 3. Monitoring
- [ ] Pipeline operations have timing metrics
- [ ] Stats/counters for processed items (samples, variants)
- [ ] Progress indicators for long operations
- [ ] Error counts and success rates tracked

### 4. No Code Duplication
- [ ] Common patterns extracted to helper functions
- [ ] Shared logic in `utils/` module
- [ ] Similar code across files consolidated

### 5. No Unused Code
- [ ] All imports are used
- [ ] No dead functions/methods
- [ ] No unused variables
- [ ] `__all__` exports match public API

---

## Phase Implementation Workflow

```
1. Review all affected files against QC checklist
2. Document findings in implementation plan
3. Fix QC issues FIRST before adding new features
4. Implement phase changes
5. Verify all tests pass
6. Re-check QC standards
7. Commit with detailed message
```

