# Nextflow Pipeline

Run gbcms at scale on HPC clusters with automatic parallelization.

## Overview

The Nextflow workflow provides:

- **Automatic parallelization** across samples
- **SLURM/HPC integration** with resource management  
- **Containerization** with Docker/Singularity
- **Resume capability** for failed runs

## Pipeline Architecture

```mermaid
flowchart TD
    CSV([📄 samplesheet.csv]):::input --> Parse[Parse Samplesheet]
    MAF([📄 variants.maf]):::input --> FilterCheck{filter_by_sample<br/>AND .maf input?}
    Parse --> FilterCheck

    FilterCheck -->|Yes| FilterMAF["FILTER_MAF<br/>(per-sample MAF extraction)"]
    FilterCheck -->|No| Ready[All samples get full variants file]

    FilterMAF --> HasData{Variants found?}
    HasData -->|Yes| Ready2[Join filtered MAF with BAM]
    HasData -->|No| Skip([⚪ Skip sample]):::skip
    FilterMAF --> Summary["PIPELINE_SUMMARY<br/>(aggregate filter stats)"]

    Ready --> Run["GBCMS_RUN<br/>(per-sample counting)"]:::run
    Ready2 --> Run

    Run --> Output([📊 VCF/MAF output]):::output
    Summary --> SummaryOut([📋 pipeline_summary.tsv]):::output

    classDef input fill:#3498db,color:#fff,stroke:#2471a3,stroke-width:2px;
    classDef run fill:#27ae60,color:#fff,stroke:#1e8449,stroke-width:2px;
    classDef output fill:#9b59b6,color:#fff,stroke:#7d3c98,stroke-width:2px;
    classDef skip fill:#95a5a6,color:#fff,stroke:#7f8c8d,stroke-width:2px;
```

## Quick Start

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile docker
```

## Documentation

| Page | Description |
|:-----|:------------|
| [Samplesheet](samplesheet.md) | Input CSV format |
| [Parameters](parameters.md) | All configuration options |
| [Examples](examples.md) | Common usage patterns |

## Related

- [CLI Reference](../cli/index.md) — For processing few samples
- [Troubleshooting](../resources/troubleshooting.md) — Common issues
