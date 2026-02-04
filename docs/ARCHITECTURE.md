# Architecture

This document describes the py-gbcms system architecture, including the Python/Rust split, data flow, and coordinate handling.

## System Overview

py-gbcms is a hybrid Python/Rust application that counts allele-supporting reads at variant positions in BAM files.

```mermaid
flowchart TB
    subgraph Python["🐍 Python Layer"]
        CLI[CLI Entry<br>cli.py] --> Pipeline[Pipeline Orchestrator<br>pipeline.py]
        Pipeline --> Reader[Input Adapters<br>io/input.py]
        Pipeline --> Writer[Output Writers<br>io/output.py]
        Reader --> Kernel[Coordinate Kernel<br>core/kernel.py]
    end
    
    subgraph Rust["🦀 Rust Layer (gbcms_rs)"]
        Counter[count_bam<br>counting.rs] --> CIGAR[CIGAR Parser]
        Counter --> Stats[Strand Bias<br>stats.rs]
    end
    
    Pipeline -->|"PyO3 binding"| Counter
    Counter -->|"BaseCounts"| Pipeline
    
    style Python fill:#3776ab,color:#fff
    style Rust fill:#dea584,color:#000
```

---

## Data Flow

```mermaid
flowchart LR
    subgraph Input
        MAF[MAF File] 
        VCF[VCF File]
        BAM[BAM Files]
        FASTA[Reference FASTA]
    end
    
    subgraph Processing
        Variants[Load Variants]
        Validate[Validate vs FASTA]
        Convert[Convert Coordinates]
        Count[Count Reads]
    end
    
    subgraph Output
        OutMAF[MAF + Counts]
        OutVCF[VCF + Counts]
    end
    
    MAF --> Variants
    VCF --> Variants
    Variants --> Validate
    FASTA --> Validate
    Validate --> Convert
    BAM --> Count
    Convert --> Count
    Count --> OutMAF
    Count --> OutVCF
```

---

## Module Organization

```
py-gbcms/
├── src/gbcms/
│   ├── cli.py           # Typer CLI entry point
│   ├── pipeline.py      # Orchestrates counting workflow
│   ├── core/
│   │   └── kernel.py    # Coordinate normalization
│   ├── io/
│   │   ├── input.py     # VCF/MAF readers
│   │   └── output.py    # VCF/MAF writers
│   ├── models/
│   │   └── core.py      # Pydantic config models
│   └── utils/
│       └── logging.py   # Structured logging
├── rust/                # Rust crate (gbcms_rs)
│   └── src/
│       ├── lib.rs       # PyO3 module entry
│       ├── counting.rs  # BAM counting logic
│       ├── stats.rs     # Fisher's exact test
│       └── types.rs     # Variant, BaseCounts
```

---

## Coordinate System

All coordinates are normalized to **0-based, half-open** internally:

```mermaid
flowchart LR
    VCF["VCF<br>(1-based)"] -->|"-1"| Internal["Internal<br>(0-based)"]
    MAF["MAF<br>(1-based)"] -->|"-1"| Internal
    Internal -->|"to Rust"| Rust["gbcms_rs"]
    Rust -->|"+1"| Output["Output<br>(1-based)"]
```

| Format | Position System | Example |
|:-------|:----------------|:--------|
| VCF | 1-based | chr1:100 |
| MAF | 1-based | chr1:100 |
| Internal | 0-based | chr1:99 |
| Rust | 0-based | chr1:99 |

---

## Rust-Python Interface

The Rust `count_bam` function is exposed via PyO3:

```python
# Python call
from gbcms_rs import count_bam, Variant

results = count_bam(
    bam_path="sample.bam",
    variants=[Variant("chr1", 99, "A", "T", "SNP")],
    fasta_path="ref.fa",
    min_mapq=20,
    min_baseq=0,
)
```

Returns `List[BaseCounts]` with:
- `ref_count`, `alt_count` (read-level)
- `ref_count_fragment`, `alt_count_fragment` (fragment-level)
- Strand-specific counts
- Fisher's exact test p-value for strand bias

---

## Configuration Hierarchy

```mermaid
flowchart TB
    GbcmsConfig --> OutputConfig
    GbcmsConfig --> ReadFilters
    GbcmsConfig --> QualityThresholds
    
    OutputConfig --> output_dir
    OutputConfig --> format
    OutputConfig --> suffix
    
    ReadFilters --> exclude_supplementary
    ReadFilters --> exclude_secondary
    ReadFilters --> count_orphans
    
    QualityThresholds --> min_mapping_quality
    QualityThresholds --> min_base_quality
```

See `src/gbcms/models/core.py` for the Pydantic model definitions.
