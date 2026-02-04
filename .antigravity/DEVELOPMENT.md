# py-gbcms Development Guide

## Quick Start

```bash
# Install in editable mode
pip install -e .

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
├── src/gbcms_rs/       # Rust extension
│   ├── src/
│   │   ├── counting.rs # BAM processing
│   │   ├── types.rs    # PyO3 bindings
│   │   └── stats.rs    # Fisher's test
│   └── Cargo.toml
├── tests/              # Test suite
└── docs/               # Documentation
```

## Building Rust Extension

```bash
cd src/gbcms_rs
maturin develop --release
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

