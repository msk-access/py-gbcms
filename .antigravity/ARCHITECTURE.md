# gbcms Architecture

## System Overview

gbcms is a **Python/Rust hybrid** tool for counting alleles at variant positions across BAM files. Python handles CLI, orchestration, and I/O; Rust handles all BAM traversal, read classification, fragment tracking, Fisher tests, mFSD analysis, and Parquet writing.

## Layer Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                     CLI Layer (Typer)                       │
│  cli.py — argument parsing, 4-layer validation, logging     │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                   Orchestration Layer                        │
│  pipeline.py — sample iteration, progress, Parquet dispatch  │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                       I/O Layer                              │
│  io/input.py  — VcfReader, MafReader, ReferenceChecker      │
│  io/output.py — VcfWriter, MafWriter  (mFSD-gated columns)  │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────┐
│                     Core / Models                            │
│  core/kernel.py — CoordinateKernel (VCF/MAF → 0-based)      │
│  models/core.py — GbcmsConfig, OutputConfig, AlignmentConfig │
└──────────────────────────┬──────────────────────────────────┘
                           │
┌─────────────────────────────────────────────────────────────┐
│                   Rust Engine (gbcms._rs)                    │
│                                                             │
│  counting/                                                   │
│    engine.rs         — Main loop, DP gating, read iteration  │
│    variant_checks.rs — check_snp/mnp/ins/del/complex         │
│    alignment.rs      — Smith-Waterman backend                │
│    pairhmm.rs        — PairHMM backend (--alignment-backend) │
│    fragment.rs       — FragmentEvidence, QNAME hashing       │
│    mfsd.rs           — Mutant Fragment Size Distribution      │
│    utils.rs          — Reconstruction, soft-clip helpers     │
│                                                             │
│  normalize/                                                  │
│    engine.rs         — Normalization pipeline                │
│    left_align.rs     — bcftools-style left-alignment         │
│    decomp.rs         — Homopolymer decomposition detection   │
│    fasta.rs          — Reference sequence fetcher            │
│    repeat.rs         — Tandem repeat detection               │
│                                                             │
│  parquet_writer.rs   — write_fsd_parquet() via ZSTD Parquet  │
│  types.rs            — Variant, BaseCounts (PyO3 bindings)   │
│  stats.rs            — Fisher's exact test                   │
└─────────────────────────────────────────────────────────────┘
```

## Data Flow

```
VCF/MAF → VcfReader/MafReader → CoordinateKernel → Variant[]
                                                       ↓
BAM + Variant[] → gbcms._rs.count_bam() → BaseCounts[]
                                              ↓
BaseCounts[] + Variant[] → VcfWriter/MafWriter → output.vcf / output.maf
                                              ↓ (--mfsd-parquet only)
BaseCounts[].ref_sizes/.alt_sizes → write_fsd_parquet() → sample.fsd.parquet
```

## Module Reference

| Module | Purpose |
|:-------|:--------|
| `cli.py` | CLI entry, 4-layer validation, logging setup |
| `pipeline.py` | Orchestration, progress, Parquet dispatch |
| `normalize.py` | Standalone normalization workflow (no counting) |
| `io/input.py` | VcfReader, MafReader, ReferenceChecker |
| `io/output.py` | VcfWriter, MafWriter with mFSD column gating |
| `core/kernel.py` | CoordinateKernel — VCF/MAF → 0-based conversion |
| `models/core.py` | GbcmsConfig, OutputConfig, AlignmentConfig (Pydantic) |
| `utils/logging.py` | Structured logging setup |
| `rust/counting/engine.rs` | Main counting loop, DP gating, Rayon parallelism |
| `rust/counting/variant_checks.rs` | check_snp/mnp/ins/del/complex, windowed scan |
| `rust/counting/alignment.rs` | Smith-Waterman implementation |
| `rust/counting/pairhmm.rs` | PairHMM backend, LLR scoring |
| `rust/counting/fragment.rs` | FragmentEvidence, quality-weighted consensus |
| `rust/counting/mfsd.rs` | Fragment size distribution analysis (KS test, LLR) |
| `rust/parquet_writer.rs` | write_fsd_parquet(), Arrow/ZSTD native Parquet |
| `rust/normalize/` | Left-alignment, decomp, fasta, repeat (7 modules) |
| `rust/types.rs` | Variant, BaseCounts PyO3 bindings |
| `rust/stats.rs` | Fisher's exact test |

## Key Design Decisions

1. **Rust for counting**: rust-htslib for BAM access; Rayon for per-variant parallelism.
2. **Python for orchestration**: Typer CLI, Pydantic config, Rich progress.
3. **0-based internal coordinates**: All coordinates converted at the Python/Rust boundary; VCF/MAF positions are 1-based externally.
4. **Reader/Writer pattern**: Pluggable adapters in `io/`. Format selection happens at writer construction, not inside the counting engine.
5. **mFSD is opt-in** (`--mfsd`): MAF/VCF writers gate 34 MAF columns and 7 VCF INFO fields behind `self.mfsd`. When flag is off, columns are absent (not NA-filled).
6. **Rust-native Parquet** (`--mfsd-parquet`): `write_fsd_parquet()` in Rust via `arrow`/`parquet` crates with ZSTD(1) compression. No `pyarrow` dependency. `ref_sizes`/`alt_sizes` are internal Rust fields — not exposed to Python via `#[pyo3(get)]`.
7. **4-layer CLI validation**: Parse-time (Typer) → Pre-model (cli.py body) → Model-time (Pydantic) → No silent skips. See `cli.py` module docstring.
8. **Fragment counting always on**: Quality-weighted consensus with QNAME hashing; discarded (ambiguous) fragments counted in DPF but not RDF/ADF.
9. **Windowed indel detection**: ±5bp scan (expanding to `max(5, repeat_span + 2)` in repeats) with 3-layer safeguards: sequence identity, closest-match tiebreak, reference context validation.
10. **Dual alignment backends**: Smith-Waterman (default) or PairHMM (`--alignment-backend hmm`). Dual-trigger local SW for borderline SW calls.

## Critical Algorithmic Features

| Feature | Location | Description |
|:--------|:---------|:------------|
| Windowed indel scan | `variant_checks.rs` | ±5bp + repeat expansion; 3-layer safeguards |
| Masked complex comparison | `variant_checks.rs` | Quality-aware masking; 3-case comparison |
| Phase 2.5 edit distance | `variant_checks.rs` | Levenshtein fallback, >1 edit margin guard |
| Phase 3 alignment | `alignment.rs`, `pairhmm.rs` | Dual-trigger SW or PairHMM |
| Fragment consensus | `fragment.rs` | u64 QNAME hash, quality-weighted, discard ambiguous |
| Multi-allelic exclusion | `engine.rs` | Sibling ALT exclusion at overlapping loci |
| Adaptive context | `normalize/repeat.rs` | Dynamic padding in tandem repeat regions |
| Strand bias | `stats.rs` | Fisher's exact test at read + fragment level |
| mFSD analysis | `counting/mfsd.rs` | KS test, LLR, pairwise comparisons across 4 fragment classes |
| Native Parquet | `parquet_writer.rs` | Arrow ListArray for variable-length size arrays, ZSTD compression |

## Configuration Model

All settings flow through `GbcmsConfig` (Pydantic):
- `OutputConfig`: `output_dir`, `format`, `suffix`, `column_prefix`, `preserve_barcode`, `show_normalization`, **`mfsd`**, **`mfsd_parquet`**
- `ReadFilters`: duplicate/secondary/supplementary/qc/improper/indel flags
- Quality thresholds: `min_mapq` (20), `min_baseq` (20), `fragment_qual_threshold` (10), `context_padding` (5)
- `AlignmentConfig`: `backend` (sw|hmm), `llr_threshold`, `gap_open_prob`, `gap_extend_prob`, repeat variants

## Test Suite (15 files)

| File | What It Tests |
|:-----|:-------------|
| `test_accuracy.py` | SNP, insertion, deletion, complex, MNP accuracy; DP invariant |
| `test_shifted_indels.py` | Windowed indel detection (15 cases) |
| `test_fuzzy_complex.py` | Quality-aware masked complex matching (14+ cases) |
| `test_filters.py` | Read filter flags |
| `test_strand_counts.py` | Strand-specific counting |
| `test_alignment_backend.py` | SW vs PairHMM concordance + integration |
| `test_fragment_consensus.py` | Fragment quality consensus, DPF invariant |
| `test_multi_allelic.py` | Sibling ALT exclusion |
| `test_dp_neither.py` | DP includes neither/third-allele reads |
| `test_normalization.py` | Left-alignment, REF validation, window expansion |
| `test_cli_sample_id.py` | CLI argument parsing, validation, error paths |
| `test_maf_reader.py` | MAF input parsing |
| `test_maf_preservation.py` | MAF column preservation |
| `test_pipeline_v2.py` | End-to-end integration |
| `test_mfsd_flag.py` | mFSD flag gating, column counts, Parquet output (8 tests) |

## Invariants (Always Assert in Tests)

- `dp >= rd + ad` — DP includes 'neither' reads
- `dpf >= rdf + adf` — DPF includes discarded ambiguous fragments
- `rd == rd_fwd + rd_rev` — strand consistency
- `ad == ad_fwd + ad_rev` — strand consistency
- MAF column count: **145** without `--mfsd`, **179** with `--mfsd` (34 additional)
- VCF INFO fields: 7 `MFSD_*` fields added only when `--mfsd` is set
