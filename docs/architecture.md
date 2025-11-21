# Architecture

## Overview

py-gbcms 2.0 uses a hybrid Python/Rust architecture for optimal performance and usability.

```
┌─────────────────────────────────────────┐
│          Python Layer (CLI)             │
│  - Typer CLI                            │
│  - Rich UI                              │
│  - Input validation                     │
│  - Output formatting                    │
└────────────┬────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────┐
│      Python I/O Layer (pysam)           │
│  - VCF/MAF reading                      │
│  - Reference validation                 │
│  - VCF/MAF writing                      │
└────────────┬────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────┐
│    Rust Counting Engine (gbcms_rs)     │
│  - BAM access (rust-htslib)             │
│  - Parallel processing (rayon)          │
│  - Fisher's exact test (statrs)         │
│  - Fragment counting                    │
└─────────────────────────────────────────┘
```

## Components

### 1. Python CLI (`src/gbcms/`)

**Modules:**
- `cli.py` - Typer-based command-line interface
- `pipeline.py` - Main orchestration logic
- `models/core.py` - Pydantic configuration models

**Responsibilities:**
- Parse command-line arguments
- Validate configuration
- Coordinate input/output
- Display progress to user

### 2. Python I/O (`src/gbcms/io/`)

**Modules:**
- `input.py` - VCF/MAF readers, reference validation
- `output.py` - VCF/MAF writers

**Key Classes:**
- `VcfReader` - Parse VCF files via pysam
- `MafReader` - Parse MAF files, normalize indels
- `ReferenceChecker` - Validate variants against FASTA
- `VcfWriter` - Write VCF with custom FORMAT fields
- `MafWriter` - Write MAF with preserved columns

### 3. Rust Counting Engine (`src/gbcms_rs/`)

**Modules:**
- `lib.rs` - PyO3 Python bindings
- `types.rs` - Variant and BaseCounts structs
- `counting.rs` - Core counting logic
- `stats.rs` - Fisher's exact test

**Key Functions:**
- `count_bam()` - Main entry point, parallelizes over variants
- `count_single_variant()` - Process one variant
- `check_allele()` - Determine if read matches ref/alt
- `fisher_exact_test()` - Calculate strand bias statistics

## Data Flow

### Input Processing

```
VCF/MAF file
    │
    ├─→ Parse with pysam
    │
    ├─→ Convert to internal Variant model (0-based)
    │
    └─→ Validate against FASTA reference
```

### Counting

```
Variants (Python list)
    │
    ├─→ Convert to Rust Variant structs
    │
    ├─→ Pass to count_bam() with BAM path
    │
    ├─→ Rayon parallel iteration over variants
    │   │
    │   ├─→ Thread 1: Process variants 0, 4, 8...
    │   ├─→ Thread 2: Process variants 1, 5, 9...
    │   ├─→ Thread 3: Process variants 2, 6, 10...
    │   └─→ Thread 4: Process variants 3, 7, 11...
    │
    └─→ Return BaseCounts list to Python
```

### Output Writing

```
BaseCounts (Rust structs)
    │
    ├─→ Convert back to Python
    │
    ├─→ Calculate VAF/FAF
    │
    ├─→ Format for VCF/MAF
    │
    └─→ Write to file
```

## Threading Model

### Data Parallelism

py-gbcms uses **data parallelism over variants**:

1. The variant list is split across threads
2. Each thread initializes its own BAM reader
3. Each thread processes its assigned variants independently
4. No shared state between threads (thread-safe)

### Why Per-Thread BAM Readers?

- **Thread Safety**: rust-htslib's `IndexedReader` is not `Sync`
- **Random Access**: Each variant may be on a different chromosome
- **Independence**: No coordination needed between threads
- **Correctness**: Avoids race conditions

### Configuration

```bash
# Default: single-threaded (safest)
--threads 1

# Multi-threaded (4 workers)
--threads 4

# Use all cores
--threads 0  # Auto-detect CPU count
```

## Memory Usage

### Per-Sample Memory

- **BAM Index**: Loaded once per thread (~10-50 MB)
- **Variant List**: Shared across threads (~1 KB per variant)
- **Read Buffer**: Per-thread, minimal (rust-htslib manages)

### Estimation

```
Memory ≈ (BAM Index Size × Threads) + (Variant Count × 1 KB)
```

**Example:**
- 10,000 variants
- 4 threads
- 30 MB BAM index

```
Memory = (30 MB × 4) + (10K × 1 KB)
       = 120 MB + 10 MB
       = 130 MB
```

## Coordinate Systems

### Internal: 0-based

All internal operations use **0-based, half-open** coordinates:
- Position 0 is the first base
- Ranges are `[start, end)` (end is exclusive)

### VCF Output: 1-based

VCF format uses **1-based, inclusive** coordinates:
- Position 1 is the first base
- Conversion: `vcf_pos = internal_pos + 1`

### MAF Output: 1-based

MAF format uses **1-based, inclusive** coordinates:
- `Start_Position` = `internal_pos + 1`
- `End_Position` = `internal_pos + len(ref)`

## Error Handling

### Python Layer

- Input validation via Pydantic
- File existence checks
- User-friendly error messages

### Rust Layer

- Result types for error propagation
- Convert Rust errors to Python exceptions
- Detailed error context

## Performance Optimizations

1. **Rust for Hot Path**: Counting loop in Rust
2. **Parallel Processing**: Rayon thread pool
3. **Efficient BAM Access**: rust-htslib with htslib C library
4. **Minimal Copies**: Zero-copy where possible
5. **Release Builds**: Optimized compilation (`--release`)

## Testing Strategy

### Unit Tests (Rust)

- Allele matching logic
- Fisher's exact test
- Coordinate conversions

### Integration Tests (Python)

- End-to-end pipeline
- VCF/MAF parsing and writing
- Filter validation
- Multi-sample processing

### Accuracy Tests

- Synthetic BAM files with known counts
- Validate exact counts for SNPs, indels, MNPs
- Filter behavior verification
