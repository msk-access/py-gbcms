# Architecture and Module Relationships

## Overview

GetBaseCounts is organized into distinct modules with clear responsibilities. This document explains how each module connects and when to use each component.

## Module Hierarchy

```
┌─────────────────────────────────────────────────────────────────┐
│                         CLI Layer                                │
│  cli.py - User interface with Typer/Rich                        │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Configuration Layer                           │
│  models.py - Pydantic models (type-safe config)                 │
│  config.py - Legacy config (backward compatibility)             │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Processing Layer                              │
│  processor.py - Main orchestration                              │
└─────┬──────────────┬──────────────┬──────────────┬─────────────┘
      │              │              │              │
      ▼              ▼              ▼              ▼
┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐
│ variant  │  │reference │  │ counter  │  │ parallel │
│  .py     │  │  .py     │  │  .py     │  │  .py     │
└──────────┘  └──────────┘  └────┬─────┘  └──────────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │numba_counter │
                           │    .py       │
                           └──────────────┘
                                  │
                                  ▼
                           ┌──────────────┐
                           │   output     │
                           │    .py       │
                           └──────────────┘
```

## Module Descriptions

### 1. CLI Layer

#### `cli.py`
**Purpose**: User interface and command-line argument parsing

**Responsibilities**:
- Parse command-line arguments with Typer
- Display beautiful output with Rich
- Validate user inputs
- Call processor with configuration

**Key Functions**:
- `count_run()` - Main counting command
- `validate_files()` - File validation
- `show_version()` - Version info
- `show_info()` - Tool information

**Dependencies**: `typer`, `rich`, `models.py`, `processor.py`

---

### 2. Configuration Layer

#### `models.py` (NEW - Pydantic)
**Purpose**: Type-safe configuration with runtime validation

**Key Classes**:
```python
GetBaseCountsConfig    # Main configuration
├── BamFileConfig      # BAM file validation
├── VariantFileConfig  # Variant file validation
├── QualityFilters     # Quality filter settings
├── OutputOptions      # Output configuration
└── PerformanceConfig  # Threading/backend settings

VariantModel           # Type-safe variant representation
VariantCounts          # Type-safe count storage
```

**When to Use**: 
- ✅ New code (recommended)
- ✅ When you need validation
- ✅ When you want type safety

#### `config.py` (LEGACY)
**Purpose**: Backward compatibility with original dataclass-based config

**When to Use**:
- ⚠️ Existing code that hasn't migrated
- ⚠️ Will be deprecated in future versions

**Migration Path**: Use `models.py` instead

---

### 3. Processing Layer

#### `processor.py`
**Purpose**: Main orchestration - coordinates all components

**Key Class**: `VariantProcessor`

**Workflow**:
```python
1. Load reference sequence (reference.py)
2. Load variants (variant.py)
3. Sort and index variants
4. For each BAM file:
   a. Create variant blocks
   b. Use parallel.py to distribute work
   c. Use counter.py to count bases
5. Write output (output.py)
```

**Dependencies**: All other modules

---

### 4. Data Loading

#### `variant.py`
**Purpose**: Load and parse variant files

**Key Classes**:
- `VariantEntry` - Legacy variant representation
- `VariantLoader` - VCF/MAF file parser

**Supports**:
- VCF format
- MAF format
- Format conversion (MAF → VCF coordinates)

**Used By**: `processor.py`

#### `reference.py`
**Purpose**: Reference sequence access

**Key Class**: `ReferenceSequence`

**Features**:
- Lazy loading with pysam.FastaFile
- Context manager support
- Base and sequence retrieval

**Used By**: `processor.py`, `variant.py` (for MAF conversion)

---

### 5. Counting Layer (THE KEY DISTINCTION)

#### `counter.py` (Pure Python Implementation)
**Purpose**: Standard base counting with pysam

**Key Class**: `BaseCounter`

**Methods**:
```python
count_bases_snp()       # Count SNP variants (DMP method)
count_bases_dnp()       # Count DNP variants  
count_bases_indel()     # Count indel variants (DMP method)
count_bases_generic()   # Generic counting for all types ⭐
filter_alignment()      # Filter reads
```

**Counting Algorithms**:

1. **DMP (Depth at Match Position) - Default**
   - Specialized methods for each variant type
   - Faster for simple variants
   - Standard algorithm from C++ version

2. **Generic Counting - Optional** (`--generic-counting`)
   - Single algorithm for all variant types
   - Parses CIGAR to extract alignment allele
   - Compares directly to ref/alt
   - Better for complex variants
   - May give slightly different results
   - Equivalent to C++ `baseCountGENERIC()`

**Characteristics**:
- ✅ Pure Python - easy to debug
- ✅ Works with pysam objects directly
- ✅ Flexible - easy to modify
- ✅ Two counting algorithms (DMP + Generic)
- ❌ Slower (baseline performance)
- ❌ No JIT compilation

**When to Use**:
- Small datasets (<10K variants)
- Debugging/development
- When Numba not available
- When you need to modify counting logic
- Complex variants (use `--generic-counting`)

**Performance**: 1x (baseline)

---

#### `numba_counter.py` (Optimized Implementation)
**Purpose**: High-performance counting with Numba JIT compilation

**Key Functions**:
```python
@jit(nopython=True, cache=True)
count_snp_base()              # Single SNP (50-100x faster)

@jit(nopython=True, parallel=True)
count_snp_batch()             # Batch SNPs (parallel)

@jit(nopython=True, cache=True)
filter_alignment_numba()      # Fast filtering

@jit(nopython=True, parallel=True)
filter_alignments_batch()     # Batch filtering

calculate_fragment_counts()   # Fragment counting
find_cigar_position()         # CIGAR parsing
```

**Characteristics**:
- ✅ 50-100x faster than pure Python
- ✅ Parallel processing with `prange`
- ✅ Cached compilation
- ❌ First call is slow (compilation)
- ❌ Requires NumPy arrays (not pysam objects)
- ❌ Harder to debug (compiled code)

**When to Use**:
- Large datasets (>10K variants)
- Production workloads
- When performance matters
- Batch processing

**Performance**: 50-100x faster than `counter.py`

---

### 6. Parallelization Layer

#### `parallel.py`
**Purpose**: Distribute work across cores/nodes

**Key Classes**:
```python
ParallelProcessor      # Unified interface
├── joblib backend     # Local parallelization
└── Ray backend        # Distributed computing

BatchProcessor         # Process in batches
VariantCounterActor    # Ray actor (stateful)
```

**Backends**:
1. **joblib** (default)
   - Local multi-core
   - Multiple backends (loky, threading, multiprocessing)
   - Best for single machine

2. **Ray** (optional)
   - Distributed across nodes
   - Fault tolerant
   - Best for clusters

**Used By**: `processor.py`

---

### 7. Output Layer

#### `output.py`
**Purpose**: Format and write results

**Key Class**: `OutputFormatter`

**Methods**:
```python
write_vcf_output()      # VCF-like format
write_maf_output()      # MAF format
write_fillout_output()  # Extended MAF with all samples
```

**Used By**: `processor.py`

---

## Data Flow

### Complete Workflow

```
User Command (cli.py)
    │
    ├─> Parse arguments
    ├─> Create GetBaseCountsConfig (models.py)
    └─> Call VariantProcessor (processor.py)
            │
            ├─> Load reference (reference.py)
            │       └─> pysam.FastaFile
            │
            ├─> Load variants (variant.py)
            │       ├─> Parse VCF/MAF
            │       └─> Create VariantEntry objects
            │
            ├─> Sort and index variants
            │
            ├─> For each BAM file:
            │   │
            │   ├─> Create variant blocks
            │   │
            │   ├─> Parallel processing (parallel.py)
            │   │       ├─> joblib: Local threads
            │   │       └─> Ray: Distributed workers
            │   │
            │   └─> For each block:
            │       │
            │       ├─> Fetch alignments (pysam)
            │       │
            │       ├─> Filter alignments
            │       │   ├─> counter.py: Python loops
            │       │   └─> numba_counter.py: JIT compiled
            │       │
            │       └─> Count bases
            │           ├─> counter.py: Pure Python
            │           │   ├─> count_bases_snp()
            │           │   ├─> count_bases_dnp()
            │           │   └─> count_bases_indel()
            │           │
            │           └─> numba_counter.py: Optimized
            │               ├─> count_snp_base()
            │               ├─> count_snp_batch()
            │               └─> filter_alignments_batch()
            │
            └─> Write output (output.py)
                    ├─> write_vcf_output()
                    ├─> write_maf_output()
                    └─> write_fillout_output()
```

---

## counter.py vs numba_counter.py

### Detailed Comparison

| Aspect | counter.py | numba_counter.py |
|--------|-----------|------------------|
| **Implementation** | Pure Python | Numba JIT compiled |
| **Input** | pysam objects | NumPy arrays |
| **Speed** | 1x (baseline) | 50-100x faster |
| **Compilation** | None | First call compiles |
| **Parallelization** | Threading only | True parallel with prange |
| **Debugging** | Easy | Harder (compiled) |
| **Flexibility** | Very flexible | Less flexible |
| **Dependencies** | pysam only | pysam + numba |
| **Use Case** | Development, small data | Production, large data |

### When Each is Used

#### `counter.py` is used when:
```python
# Default for small datasets
if num_variants < 10000 and not config.use_numba:
    counter = BaseCounter(config)
    counts = counter.count_bases_snp(variant, alignments, sample)
```

#### `numba_counter.py` is used when:
```python
# Enabled by default for performance
if config.use_numba:  # Default: True
    # Convert pysam data to NumPy arrays
    bases_array = np.array([aln.query_sequence for aln in alignments])
    quals_array = np.array([aln.query_qualities for aln in alignments])
    
    # Use Numba-optimized counting
    from numba_counter import count_snp_batch
    counts = count_snp_batch(bases_array, quals_array, ...)
```

### Integration Pattern

The `processor.py` can use both:

```python
class VariantProcessor:
    def __init__(self, config):
        self.config = config
        self.counter = BaseCounter(config)  # Always available
        self.use_numba = config.performance.use_numba
    
    def count_variant_block(self, variants, alignments):
        if self.use_numba and len(variants) > 100:
            # Use Numba for large batches
            return self._count_with_numba(variants, alignments)
        else:
            # Use pure Python for small batches
            return self._count_with_python(variants, alignments)
    
    def _count_with_python(self, variants, alignments):
        """Use counter.py"""
        for variant in variants:
            self.counter.count_bases_snp(variant, alignments, sample)
    
    def _count_with_numba(self, variants, alignments):
        """Use numba_counter.py"""
        # Convert to NumPy arrays
        data = self._prepare_numba_data(alignments)
        # Batch process with Numba
        from numba_counter import count_snp_batch
        return count_snp_batch(**data)
```

---

## Configuration Flow

### Using Pydantic Models (Recommended)

```python
from getbasecounts.models import GetBaseCountsConfig, PerformanceConfig

config = GetBaseCountsConfig(
    fasta_file=Path("ref.fa"),
    bam_files=[...],
    variant_files=[...],
    performance=PerformanceConfig(
        use_numba=True,        # Use numba_counter.py
        num_threads=16,
        backend='joblib',
    ),
)

processor = VariantProcessor(config)
processor.process()
```

### Legacy Config (Backward Compatible)

```python
from getbasecounts.config import Config

config = Config(
    fasta_file="ref.fa",
    bam_files={...},
    variant_files=[...],
    num_threads=16,
)

processor = VariantProcessor(config)
processor.process()
```

---

## Performance Optimization Path

### Level 1: Pure Python (counter.py)
```python
config = GetBaseCountsConfig(
    performance=PerformanceConfig(
        use_numba=False,
        num_threads=1,
    )
)
# Speed: 1x
```

### Level 2: Multi-threaded Python (counter.py + joblib)
```python
config = GetBaseCountsConfig(
    performance=PerformanceConfig(
        use_numba=False,
        num_threads=8,
        backend='joblib',
    )
)
# Speed: ~4-6x
```

### Level 3: Numba Single-threaded (numba_counter.py)
```python
config = GetBaseCountsConfig(
    performance=PerformanceConfig(
        use_numba=True,
        num_threads=1,
    )
)
# Speed: ~50-100x
```

### Level 4: Numba + joblib (numba_counter.py + parallel.py)
```python
config = GetBaseCountsConfig(
    performance=PerformanceConfig(
        use_numba=True,
        num_threads=16,
        backend='joblib',
    )
)
# Speed: ~200-400x
```

### Level 5: Numba + Ray (numba_counter.py + Ray)
```python
config = GetBaseCountsConfig(
    performance=PerformanceConfig(
        use_numba=True,
        num_threads=32,
        backend='ray',
        use_ray=True,
    )
)
# Speed: ~500-1000x (on cluster)
```

---

## Testing Strategy

### Unit Tests

Each module has its own tests:
- `tests/test_config.py` - Configuration
- `tests/test_variant.py` - Variant loading
- `tests/test_counter.py` - **counter.py** functions
- `tests/test_reference.py` - Reference access
- `tests/test_output.py` - Output formatting
- `tests/test_cli.py` - CLI interface

### Integration Tests

- `scripts/test_vcf_workflow.sh` - End-to-end VCF
- `scripts/test_maf_workflow.sh` - End-to-end MAF

---

## Summary

### Key Takeaways

1. **counter.py** = Pure Python, flexible, slower
2. **numba_counter.py** = JIT compiled, fast, less flexible
3. **Both can coexist** - use based on workload
4. **processor.py** orchestrates everything
5. **models.py** provides type safety
6. **parallel.py** handles distribution

### Decision Tree

```
Need to count bases?
├─> Small dataset (<10K variants)
│   └─> Use counter.py
│
├─> Large dataset (>10K variants)
│   └─> Use numba_counter.py
│
├─> Need to debug/modify counting logic
│   └─> Use counter.py
│
└─> Production workload
    └─> Use numba_counter.py + parallel.py
```

### Module Dependencies

```
cli.py
 └─> models.py
      └─> processor.py
           ├─> variant.py
           ├─> reference.py
           ├─> counter.py (optional)
           ├─> numba_counter.py (optional)
           ├─> parallel.py
           └─> output.py
```

All modules are designed to work together seamlessly while maintaining clear separation of concerns! 🎯
