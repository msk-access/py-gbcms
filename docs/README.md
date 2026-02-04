# py-gbcms

> **Get Base Counts Multi-Sample** — High-performance variant counting from BAM files

[![Version](https://img.shields.io/pypi/v/py-gbcms)](https://pypi.org/project/py-gbcms/)
[![Python](https://img.shields.io/pypi/pyversions/py-gbcms)](https://pypi.org/project/py-gbcms/)
[![License](https://img.shields.io/github/license/msk-access/py-gbcms)](https://github.com/msk-access/py-gbcms/blob/main/LICENSE)

## What It Does

gbcms extracts **allele counts** and **variant metrics** at specified positions in BAM files:

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

**→ [Full Installation Guide](INSTALLATION.md)** | **→ [CLI Examples](quick-start.md)**

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
| **CLI** | 1-10 samples, local/single server | [Quick Start](quick-start.md) |
| **Nextflow** | 10+ samples, HPC/SLURM | [Nextflow Guide](NEXTFLOW.md) |

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

**→ [Technical Details](ARCHITECTURE.md)**

---

## Links

- **[Installation](INSTALLATION.md)** — PyPI, Docker, source
- **[CLI Guide](quick-start.md)** — Command examples  
- **[Nextflow](NEXTFLOW.md)** — HPC pipeline
- **[Contributing](CONTRIBUTING.md)** — Development guide
- **[Changelog](CHANGELOG.md)** — Version history
