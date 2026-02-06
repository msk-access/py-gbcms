# py-gbcms

> **Get Base Counts Multi-Sample** — High-performance variant counting from BAM files

[![Version](https://img.shields.io/pypi/v/py-gbcms)](https://pypi.org/project/py-gbcms/)
[![Python](https://img.shields.io/pypi/pyversions/py-gbcms)](https://pypi.org/project/py-gbcms/)
[![License](https://img.shields.io/github/license/msk-access/py-gbcms)](https://github.com/msk-access/py-gbcms/blob/main/LICENSE)

## What It Does

GBCMS extracts **allele counts** and **variant metrics** at specified positions in BAM files:

```mermaid
flowchart LR
    subgraph Input
        VCF[VCF/MAF]
        BAM[BAM Files]
    end
    
    subgraph Output
        Counts[Allele Counts]
        Metrics[VAF, Strand Bias]
    end
    
    VCF --> Engine[gbcms]
    BAM --> Engine
    Engine --> Counts
    Engine --> Metrics
```

### Key Metrics

| Metric | Formula | Description |
|:-------|:--------|:------------|
| **VAF** | `AD / (RD + AD)` | Variant Allele Frequency |
| **Strand Bias** | Fisher's exact test | Detect sequencing artifacts |
| **Fragment Counts** | Deduplicated pairs | PCR-aware counting |

---

## Quick Start

```bash
# Install
pip install py-gbcms

# Run
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/
```

**→ [Full Installation Guide](getting-started/installation.md)** | **→ [CLI Examples](getting-started/quickstart.md)**

---

## Choose Your Workflow

```mermaid
flowchart TD
    Start[How many samples?] --> Few{1-10 samples}
    Start --> Many{10+ samples}
    
    Few --> CLI[Use CLI<br/>gbcms run ...]
    Many --> HPC{HPC cluster?}
    
    HPC --> |Yes| Nextflow[Use Nextflow<br/>nextflow run ...]
    HPC --> |No| CLI
    
    style CLI fill:#4caf50,color:white
    style Nextflow fill:#2196f3,color:white
```

| Workflow | Best For | Guide |
|:---------|:---------|:------|
| **CLI** | 1-10 samples, local/single server | [Quick Start](getting-started/quickstart.md) |
| **Nextflow** | 10+ samples, HPC/SLURM | [Nextflow Guide](nextflow/index.md) |

---

## Architecture

Python/Rust hybrid for maximum performance:

```mermaid
flowchart TB
    subgraph Python["🐍 Python"]
        CLI[CLI] --> Pipeline[Orchestration]
        Pipeline --> IO[VCF/MAF I/O]
    end
    
    subgraph Rust["🦀 Rust"]
        Counter[BAM Counting]
        Stats[Fisher's Test]
    end
    
    Pipeline --> Counter
    Counter --> Stats
    
    style Python fill:#3776ab,color:#fff
    style Rust fill:#dea584,color:#000
```

**→ [Technical Details](reference/architecture.md)**

---

## Documentation

| Section | Description |
|:--------|:------------|
| [Getting Started](getting-started/index.md) | Installation and first run |
| [CLI Reference](cli/index.md) | Command-line usage |
| [Nextflow Pipeline](nextflow/index.md) | HPC workflow |
| [Reference](reference/architecture.md) | Architecture and formats |
| [Development](development/developer-guide.md) | Contributing guide |
