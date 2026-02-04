# py-gbcms Architecture

## Overview

py-gbcms is a Python/Rust hybrid tool for counting bases at variant positions in BAM files.

## Layer Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     CLI Layer (Typer)                       │
│  cli.py - argument parsing, logging setup, error handling   │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                   Orchestration Layer                        │
│  pipeline.py - workflow, progress, sample iteration         │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                     I/O Layer                                │
│  io/input.py  - VcfReader, MafReader, ReferenceChecker      │
│  io/output.py - VcfWriter, MafWriter                        │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                    Core Layer                                │
│  core/kernel.py - CoordinateKernel (VCF/MAF normalization)  │
│  models/core.py - Variant, GbcmsConfig (Pydantic)           │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                   Rust Engine (gbcms._rs)                    │
│  counting.rs - BAM processing, strand bias, Fisher's test   │
│  types.rs    - Variant, BaseCounts (PyO3 bindings)          │
└─────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

1. **Rust for BAM Processing**: High-performance counting with rust-htslib
2. **Python for Orchestration**: Typer CLI, Pydantic config, rich progress
3. **Coordinate Kernel**: Centralized VCF/MAF → 0-based coordinate conversion
4. **Reader/Writer Pattern**: Pluggable I/O adapters for format handling

## Module Breakdown

| Module | LOC | Purpose |
|:-------|:----|:--------|
| `cli.py` | ~170 | CLI entry, argument parsing |
| `pipeline.py` | ~250 | Workflow orchestration |
| `io/input.py` | ~230 | Variant file readers |
| `io/output.py` | ~360 | Result file writers |
| `core/kernel.py` | ~130 | Coordinate transforms |
| `models/core.py` | ~130 | Data models |
| `utils/logging.py` | ~110 | Logging utilities |
| `rust/` | ~500 | Rust extension (bundled as gbcms._rs) |

## Data Flow

```
VCF/MAF → VcfReader/MafReader → CoordinateKernel → Variant[]
                                                       ↓
BAM + Variant[] → gbcms._rs.count_bam() → BaseCounts[]
                                              ↓
BaseCounts[] + Variant[] → VcfWriter/MafWriter → Output File
```

## Configuration

All settings flow through `GbcmsConfig` (Pydantic model):
- Input paths (variants, BAMs, reference)
- Output settings (dir, format, suffix)
- Filters (mapq, baseq, duplicates, etc.)
- Performance (threads)

## Testing

Tests in `tests/`:
- `test_accuracy.py` - Synthetic BAM accuracy tests
- `test_filters.py` - Read filter flag tests
- `test_pipeline_v2.py` - Integration tests
