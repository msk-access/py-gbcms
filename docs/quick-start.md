# Quick Start Guide

This guide shows you how to get started with `py-gbcms` using the **standalone CLI** for processing one or a few samples.

> **Processing many samples?** Use the [Nextflow Workflow](NEXTFLOW.md) instead for automatic parallelization on HPC clusters.

## Prerequisites

- Python >= 3.10
- Rust toolchain (for installation from source)
- BAM files with index (.bai)
- Reference FASTA with index (.fai)
- Variants file (VCF or MAF)

Install via `pip install py-gbcms` or see the project README for detailed setup instructions.

## Basic Usage

### Single Sample

Count variants for one sample:

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample1.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/sample1.vcf`

### Multiple Samples

Process multiple samples sequentially:

```bash
# Sample 1
gbcms run --variants variants.vcf --bam sample1.bam --fasta ref.fa --output-dir results/

# Sample 2
gbcms run --variants variants.vcf --bam sample2.bam --fasta ref.fa --output-dir results/

# Sample 3
gbcms run --variants variants.vcf --bam sample3.bam --fasta ref.fa --output-dir results/
```

Or use a BAM list file:

```bash
# Create bam_list.txt with:
# sample1 /path/to/sample1.bam
# sample2 /path/to/sample2.bam

gbcms run \
    --variants variants.vcf \
    --bam-list bam_list.txt \
    --fasta reference.fa \
    --output-dir results/
```

> **Note:** The CLI processes samples sequentially. For parallel processing of many samples, use the [Nextflow Workflow](NEXTFLOW.md).

## Common Options

### Output Format

**VCF (default):**
```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/ --format vcf
```

**MAF:**
```bash
gbcms run --variants variants.maf --bam sample.bam --fasta ref.fa --output-dir results/ --format maf
```

### Custom Sample IDs

Override the sample name:

```bash
gbcms run \
    --variants variants.vcf \
    --bam MySampleID:sample.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/MySampleID.vcf`

### Output Suffix

Add suffix to output filenames:

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/ \
    --suffix .genotyped
```

**Output:** `results/sample.genotyped.vcf`

### Threading

Use multiple threads for processing:

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/ \
    --threads 4
```

### Quality Filters

**Minimum mapping quality:**
```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/ --min-mapq 30
```

**Minimum base quality:**
```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/ --min-baseq 20
```

**Filter duplicates (default: enabled):**
```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/ --filter-duplicates
```

**Filter secondary alignments:**
```bash
gbcms run --variants variants.vcf --bam sample.bam --fasta ref.fa --output-dir results/ --filter-secondary
```

## Complete Example

Process a sample with strict filtering:

```bash
gbcms run \
    --variants variants.vcf \
    --bam TumorSample:tumor.bam \
    --fasta hg19.fa \
    --output-dir genotyped_results/ \
    --format vcf \
    --suffix .genotyped \
    --threads 8 \
    --min-mapq 30 \
    --min-baseq 20 \
    --filter-duplicates \
    --filter-secondary \
    --filter-supplementary
```

**Output:** `genotyped_results/TumorSample.genotyped.vcf`

## Using Docker

Run via Docker container:

```bash
docker run --rm -v $(pwd):/data ghcr.io/msk-access/py-gbcms:2.0.0 \
    gbcms run \
    --variants /data/variants.vcf \
    --bam /data/sample.bam \
    --fasta /data/reference.fa \
    --output-dir /data/results/
```

## Next Steps

- **Many samples on HPC:** See [Nextflow Workflow](NEXTFLOW.md)
- **Usage patterns:** See [Usage Overview](WORKFLOWS.md)
