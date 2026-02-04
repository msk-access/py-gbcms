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
