# py-gbcms Nextflow Workflow

This directory contains a Nextflow workflow for running `py-gbcms` on multiple samples in parallel.

## Prerequisites

1. [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (>=21.10.3)
2. [Docker](https://docs.docker.com/engine/installation/) or [Singularity](https://www.sylabs.io/guides/3.0/user-guide/)

## Quick Start

### 1. Prepare a samplesheet (CSV format)

**Basic samplesheet:**
```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,
```

**With per-sample suffix (for multiple BAM types per sample):**
```csv
sample,bam,bai,suffix
sample1,/path/to/sample1.duplex.bam,,-duplex
sample1,/path/to/sample1.simplex.bam,,-simplex
sample1,/path/to/sample1.unfiltered.bam,,-unfiltered
sample2,/path/to/sample2.bam,,
```

**Output:**
- `sample1-duplex.vcf`
- `sample1-simplex.vcf`  
- `sample1-unfiltered.vcf`
- `sample2.vcf` (or `sample2{--suffix}.vcf` if global suffix set)

**Notes:**
- `bai` column is optional - will auto-discover `<bam>.bai` if not provided
- `suffix` column is optional - per-row suffix overrides global `--suffix` parameter

### 2. Run the workflow

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
| `--fasta` | Path to reference FASTA (with .fai index) |

### Optional
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `results` |
| `--format` | Output format (`vcf` or `maf`) | `vcf` |
| `--suffix` | Suffix to append to output filenames | `''` (empty) |
| `--min_mapq` | Minimum mapping quality | `20` |
| `--min_baseq` | Minimum base quality | `0` |
| `--filter_duplicates` | Filter duplicate reads | `true` |
| `--filter_secondary` | Filter secondary alignments | `false` |
| `--filter_supplementary` | Filter supplementary alignments | `false` |
| `--filter_qc_failed` | Filter QC failed reads | `false` |
| `--filter_improper_pair` | Filter improperly paired reads | `false` |
| `--filter_indel` | Filter reads with indels | `false` |

### Resource Limits
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max_cpus` | Maximum CPUs per process | `16` |
| `--max_memory` | Maximum memory per process | `128.GB` |
| `--max_time` | Maximum time per process | `240.h` |

## Profiles

- `-profile docker`: Use Docker containers (recommended for local)
- `-profile singularity`: Use Singularity images (recommended for HPC)
- `-profile slurm`: Run on SLURM cluster with Singularity (queue: `cpu_medium`)
- `-profile debug`: Print hostname for debugging

## Customizing for Your Cluster

Edit `nextflow/nextflow.config` to customize the SLURM profile:

```groovy
slurm {
    process.executor       = 'slurm'
    process.queue          = 'your_queue_name'  // Change this
    singularity.enabled    = true
    singularity.autoMounts = true
    docker.enabled         = false
}
```

## Output

Results are published to `${params.outdir}/gbcms/`:
- VCF files: `<sample>.gbcms.vcf`
- MAF files: `<sample>.gbcms.maf`

Pipeline info and logs are in `${params.outdir}/pipeline_info/`.
