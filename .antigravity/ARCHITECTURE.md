# gbcms Architecture

## Overview

gbcms is a Python/Rust hybrid tool for counting bases at variant positions in BAM files.

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
┌─────────────────────────────────────────────────────────────┐
│                   Rust Engine (gbcms._rs)                    │
│  counting/                                                   │
│    engine.rs     - Main loop, DP gating, read iteration       │
│    variant_checks.rs - check_snp/mnp/ins/del/complex          │
│    alignment.rs  - Smith-Waterman implementation               │
│    pairhmm.rs    - PairHMM alignment backend                   │
│    fragment.rs   - FragmentEvidence, quality consensus          │
│    utils.rs      - Reconstruction, soft-clip helpers           │
│  normalize/                                                   │
│    engine.rs     - Normalization pipeline                      │
│    left_align.rs - bcftools-style left-alignment               │
│    decomp.rs     - Homopolymer decomposition                   │
│    fasta.rs      - Reference sequence fetcher                  │
│    repeat.rs     - Tandem repeat detection                     │
│  types.rs    - Variant, BaseCounts (PyO3 bindings)            │
│  stats.rs    - Fisher's exact test                            │
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
| `cli.py` | ~350 | CLI entry, argument parsing, AlignmentConfig |
| `pipeline.py` | ~450 | Workflow orchestration, progress |
| `io/input.py` | ~140 | Variant file readers |
| `io/output.py` | ~430 | Result file writers |
| `core/kernel.py` | ~230 | Coordinate transforms |
| `models/core.py` | ~260 | Data models (GbcmsConfig, AlignmentConfig) |
| `utils/logging.py` | ~140 | Logging utilities |
| `rust/counting/` | ~3400 | Counting engine (7 modules) |
| `rust/normalize/` | ~1400 | Normalization engine (7 modules) |
| `rust/` (other) | ~200 | lib.rs, stats.rs, types.rs |
| **Total** | **~7000** | **Python ~2100 + Rust ~5000** |

## Data Flow

```
VCF/MAF → VcfReader/MafReader → CoordinateKernel → Variant[]
                                                       ↓
BAM + Variant[] → gbcms._rs.count_bam() → BaseCounts[]
                                              ↓
BaseCounts[] + Variant[] → VcfWriter/MafWriter → Output File
```

## Key Algorithmic Features

| Feature | Module | Description |
|:--------|:-------|:------------|
| Windowed Indel Detection | `counting/variant_checks.rs` | ±5bp scan with 3-layer safeguards (sequence identity, closest match, reference context) |
| Masked Complex Comparison | `counting/variant_checks.rs` | Quality-aware masking with ambiguity detection; 3-case comparison (equal-length, ALT-only, REF-only) |
| Phase 2.5 Edit Distance | `counting/variant_checks.rs` | Levenshtein distance fallback with >1 edit margin guard; skips when ref_len > 2×alt_len |
| Phase 3 Alignment | `counting/alignment.rs`, `counting/pairhmm.rs` | Dual backend: SW (default) or PairHMM. Dual-trigger local SW for borderline calls |
| Fragment Consensus | `counting/fragment.rs` | `FragmentEvidence` struct with u64 QNAME hashing, quality-weighted consensus, discard on ambiguity |
| Multi-Allelic Exclusion | `counting/engine.rs` | Sibling ALT exclusion at overlapping loci to prevent false REF inflation |
| Adaptive Context | `normalize/repeat.rs` | Dynamic context padding expansion in tandem repeat regions |
| Strand Bias | `stats.rs` | Fisher's exact test at both read and fragment level |

## Configuration

All settings flow through `GbcmsConfig` (Pydantic model):
- Input paths (variants, BAMs, reference)
- Output settings (dir, format, suffix)
- Filters (mapq, baseq=20 default, duplicates, etc.)
- Alignment backend (sw/hmm, LLR threshold, gap probabilities)
- Performance (threads)

## Testing

14 test files across `tests/`:
- `test_accuracy.py` - Synthetic BAM accuracy (SNP, indel, complex, MNP)
- `test_shifted_indels.py` - Windowed indel detection (15 cases)
- `test_fuzzy_complex.py` - Quality-aware masked comparison (14 cases)
- `test_filters.py` - Read filter flag tests
- `test_strand_counts.py` - Strand-specific counting
- `test_alignment_backend.py` - SW vs PairHMM concordance + integration
- `test_fragment_consensus.py` - Fragment quality consensus
- `test_multi_allelic.py` - Multi-allelic exclusion
- `test_dp_neither.py` - DP/neither classification edge cases
- `test_normalization.py` - Variant normalization
- `test_cli_sample_id.py` - CLI argument parsing
- `test_maf_reader.py` / `test_maf_preservation.py` - MAF I/O
- `test_pipeline_v2.py` - Integration tests
