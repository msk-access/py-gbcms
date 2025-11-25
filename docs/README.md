# Introduction

**gbcms** (Get Base Counts Multi-Sample) is a high-performance tool for extracting base counts and variant metrics from BAM files.

## Key Features

- **Rust-Powered Engine**: Maximum performance and memory efficiency
- **Accurate Variant Handling**: Full VCF/MAF support with rigorous normalization
- **Comprehensive Metrics**: Depth, strand counts, fragment counts, VAF, and Fisher's strand bias
- **Modern CLI**: User-friendly interface with rich output
- **Flexible Deployment**: Standalone CLI or Nextflow workflow

## Two Ways to Use gbcms

### 🔧 Standalone CLI
For processing **1-10 samples** locally or on a single server.

```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/
```

**→ [CLI Quick Start](quick-start.md)**

### 🔄 Nextflow Workflow
For processing **10+ samples** in parallel on HPC clusters.

```bash
nextflow run nextflow/main.nf --input samples.csv --variants variants.vcf --fasta ref.fa -profile slurm
```

**→ [Nextflow Workflow Guide](NEXTFLOW.md)**

---

## Architecture

gbcms uses a hybrid Python/Rust architecture:
- **Python**: CLI, input parsing (VCF/MAF), orchestration, and output formatting
- **Rust**: Computationally intensive read iteration, base pileup, and statistics

---

## Next Steps

Choose your path:

**New to gbcms?**
1. [Usage Overview](WORKFLOWS.md) - Choose CLI or Nextflow
2. [CLI Quick Start](quick-start.md) OR [Nextflow Guide](NEXTFLOW.md)

**Contributing:**
- [Contributing Guide](CONTRIBUTING.md)
