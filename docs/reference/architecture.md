# Architecture

py-gbcms uses a hybrid Python/Rust architecture for maximum performance.

## System Overview

```mermaid
flowchart TB
    subgraph Python["🐍 Python Layer"]
        CLI[CLI<br/>cli.py] --> Pipeline[Orchestration<br/>pipeline.py]
        Pipeline --> Reader[Input Adapters<br/>VcfReader, MafReader]
        Pipeline --> Writer[Output Writers<br/>VcfWriter, MafWriter]
    end
    
    subgraph Rust["🦀 Rust Layer (gbcms._rs)"]
        Counter[count_bam<br/>counting.rs] --> CIGAR[CIGAR Parser]
        Counter --> Stats[Strand Bias<br/>stats.rs]
    end
    
    Pipeline -->|"PyO3"| Counter
    Counter -->|"BaseCounts"| Pipeline
    
    style Python fill:#3776ab,color:#fff
    style Rust fill:#dea584,color:#000
```

---

## Data Flow

```mermaid
flowchart LR
    subgraph Input
        VCF[VCF/MAF]
        BAM[BAM Files]
        FASTA[Reference]
    end
    
    subgraph Process
        Load[Load Variants]
        Validate[Validate vs Ref]
        Count[Count Reads]
    end
    
    subgraph Output
        Result[VCF/MAF + Counts]
    end
    
    VCF --> Load --> Validate
    FASTA --> Validate
    Validate --> Count
    BAM --> Count
    Count --> Result
```

---

## Coordinate System

All coordinates normalized to **0-based, half-open** internally:

```mermaid
flowchart LR
    VCF["VCF (1-based)"] -->|"-1"| Internal["Internal (0-based)"]
    MAF["MAF (1-based)"] -->|"-1"| Internal
    Internal -->|"to Rust"| Rust["gbcms._rs"]
    Rust -->|"+1"| Output["Output (1-based)"]
```

| Format | System | Example |
|:-------|:-------|:--------|
| VCF input | 1-based | chr1:100 |
| Internal | 0-based | chr1:99 |
| Output | 1-based | chr1:100 |

---

## Formulas

### Variant Allele Frequency (VAF)

```
VAF = AD / (RD + AD)
```

Where:
- **AD** = Alternate allele read count
- **RD** = Reference allele read count

### Strand Bias (Fisher's Exact Test)

```
         |  Forward  Reverse  |
    -----+--------------------+
    Ref  |    a        b      |
    Alt  |    c        d      |
    -----+--------------------+
    
    p-value = Fisher's exact test on 2×2 contingency table
```

Low p-value (< 0.05) indicates potential strand bias artifact.

---

## Module Structure

```
src/gbcms/
├── cli.py           # Typer CLI
├── pipeline.py      # Orchestration
├── core/
│   └── kernel.py    # Coordinate normalization
├── io/
│   ├── input.py     # VcfReader, MafReader
│   └── output.py    # VcfWriter, MafWriter
├── models/
│   └── core.py      # Pydantic config
└── utils/
    └── logging.py   # Structured logging

rust/src/
├── lib.rs           # PyO3 module (_rs)
├── counting.rs      # BAM processing
├── stats.rs         # Fisher's exact test
└── types.rs         # Variant, BaseCounts
```

---

## Configuration

All settings via `GbcmsConfig` (Pydantic model):

```mermaid
flowchart TB
    GbcmsConfig --> OutputConfig[Output Settings]
    GbcmsConfig --> ReadFilters[Read Filters]
    GbcmsConfig --> QualityThresholds[Quality Thresholds]
    
    OutputConfig --> D1[output_dir, format, suffix]
    ReadFilters --> D2[exclude_secondary, exclude_duplicates]
    QualityThresholds --> D3[min_mapq, min_baseq]
```

See [models/core.py](file:///src/gbcms/models/core.py) for definitions.

## Related

- [Variant Counting](variant-counting.md) — How each variant type is counted
- [Input Formats](input-formats.md) — VCF and MAF specifications
- [Glossary](glossary.md) — Term definitions
