# Nextflow Workflow Guide

This guide covers running `py-gbcms` as a Nextflow workflow for processing multiple samples in parallel, particularly on HPC clusters.

## Overview

The Nextflow workflow provides:
- **Automatic parallelization** across samples
- **SLURM/HPC integration** with resource management
- **Containerization** with Docker/Singularity
- **Resume capability** for failed runs
- **Reproducible pipelines**

## Prerequisites

1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) >= 21.10.3
2. One of:
   - [Docker](https://docs.docker.com/engine/installation/) (for local)
   - [Singularity](https://sylabs.io/guides/3.0/user-guide/) (for HPC)

**Install Nextflow:**
```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or any directory in your PATH
```

## Quick Start

### 1. Prepare Samplesheet

Create a CSV file with your samples:

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,
sample3,/path/to/sample3.bam,/path/to/sample3.bam.bai
```

**Or with per-sample suffix (for multiple BAM types):**
```csv
sample,bam,bai,suffix
sample1,/path/to/sample1.duplex.bam,,-duplex
sample1,/path/to/sample1.simplex.bam,,-simplex
sample1,/path/to/sample1.unfiltered.bam,,-unfiltered
sample2,/path/to/sample2.bam,,
```

**Notes:**
- `bai` column is optional - will auto-discover `<bam>.bai` if not provided
- `suffix` column is optional - per-row suffix overrides global `--suffix` parameter
- BAI files must exist or workflow will fail early with clear error

### 2. Run the Workflow

**Local with Docker:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile docker
```

**SLURM cluster with Singularity:**
```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile slurm
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV |
| `--variants` | Path to VCF/MAF variants file |
| `--fasta` | Reference FASTA (with .fai index) |

### Output Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |

### Filtering Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_mapq` | `20` | Minimum mapping quality |
| `--min_baseq` | `0` | Minimum base quality |
| `--filter_duplicates` | `true` | Filter duplicate reads |
| `--filter_secondary` | `false` | Filter secondary alignments |
| `--filter_supplementary` | `false` | Filter supplementary alignments |
| `--filter_qc_failed` | `false` | Filter QC failed reads |
| `--filter_improper_pair` | `false` | Filter improperly paired reads |
| `--filter_indel` | `false` | Filter reads with indels |

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Maximum CPUs per job |
| `--max_memory` | `128.GB` | Maximum memory per job |
| `--max_time` | `240.h` | Maximum runtime per job |

## Execution Profiles

### Docker (Local)
```bash
-profile docker
```
- Uses Docker containers
- Best for local development
- Requires Docker installed

### Singularity (HPC)
```bash
-profile singularity
```
- Uses Singularity images
- Best for HPC without SLURM
- Requires Singularity installed

### SLURM (HPC Cluster)
```bash
-profile slurm
```
- Submits jobs to SLURM
- Uses Singularity containers
- Queue: `cmobic_cpu` (customizable)

## Customizing for Your Cluster

Edit `nextflow/nextflow.config` to customize the SLURM profile:

```groovy
slurm {
    process.executor       = 'slurm'
    process.queue          = 'your_queue_name'  // Change this
    process.clusterOptions = '--account=your_account'  // Add if needed
    singularity.enabled    = true
    singularity.autoMounts = true
}
```

Common customizations:
```groovy
process {
    withName: GBCMS_RUN {
        cpus   = 8          // CPUs per sample
        memory = 16.GB      // Memory per sample
        time   = 6.h        // Time limit per sample
    }
}
```

## Output Structure

Results are organized in `${outdir}/`:

```
results/
├── gbcms/
│   ├── sample1.vcf        # Or .maf
│   ├── sample2.vcf
│   └── sample3.vcf
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## Advanced Usage

### Resume Failed Runs

Nextflow caches completed tasks. Resume from where it failed:

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile slurm \
    -resume
```

### Custom Suffix

Add suffix to output filenames:

```bash
--suffix .genotyped
# Output: sample1.genotyped.vcf
```

### MAF Output

Generate MAF instead of VCF:

```bash
--format maf
# Output: sample1.maf
```

### Strict Filtering

Enable all filters for high-quality genotyping:

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --filter_duplicates true \
    --filter_secondary true \
    --filter_supplementary true \
    --filter_qc_failed true \
    -profile slurm
```

## Monitoring

### View Running Jobs

```bash
# SLURM
squeue -u $USER

# Nextflow
nextflow log
```

### Check Progress

Nextflow prints real-time progress:
```
[c3/a1b2c3] GBCMS_RUN (sample1) [100%] 10 of 10 ✔
```

### Execution Report

After completion, view the HTML report:
```bash
open results/pipeline_info/execution_report.html
```

## Troubleshooting

### Job Failed with Error

Check the work directory in error message:
```bash
cat work/c3/a1b2c3/.command.log
```

### Out of Memory

Increase memory in config:
```groovy
process {
    withName: GBCMS_RUN {
        memory = 32.GB
    }
}
```

### Wrong Queue

Update queue name in `nextflow/nextflow.config`:
```groovy
process.queue = 'your_queue_name'
```

### Missing Container

Pull the container manually:
```bash
# Singularity
singularity pull docker://ghcr.io/msk-access/py-gbcms:2.0.0

# Docker
docker pull ghcr.io/msk-access/py-gbcms:2.0.0
```

## Comparison with CLI

| Feature | CLI | Nextflow |
|---------|-----|----------|
| Multiple samples | Sequential | Parallel |
| Resource management | Manual | Automatic |
| Retry failed jobs | Manual | Automatic |
| HPC integration | Manual scripts | Built-in |
| Resume capability | No | Yes |

**When to use CLI instead:** See [Usage Patterns](WORKFLOWS.md)

## Next Steps

- Read the [CLI Reference](CLI_FEATURES.md) for parameter details
- Check [Input & Output Formats](INPUT_OUTPUT.md) for file specifications
- See `nextflow/README.md` for additional workflow documentation
