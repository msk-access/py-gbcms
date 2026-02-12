# Nextflow Parameters

Complete reference for all pipeline parameters.

## Required Parameters

| Parameter | Description |
|:----------|:------------|
| `--input` | Path to [samplesheet CSV](samplesheet.md) |
| `--variants` | Path to [VCF/MAF](../reference/input-formats.md) variants file |
| `--fasta` | Reference FASTA (with .fai index) |

## Output Options

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--outdir` | `results` | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |
| `--column_prefix` | `''` | Prefix for gbcms count columns in MAF output |
| `--preserve_barcode` | `false` | Keep original Tumor_Sample_Barcode from input MAF |

## Filtering Options

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--min_mapq` | `20` | Minimum MAPQ |
| `--min_baseq` | `20` | Minimum BASEQ |
| `--filter_duplicates` | `true` | Filter duplicate reads |
| `--filter_secondary` | `false` | Filter secondary alignments |
| `--filter_supplementary` | `false` | Filter supplementary alignments |
| `--filter_qc_failed` | `false` | Filter QC failed reads |
| `--filter_improper_pair` | `false` | Filter improperly paired reads |
| `--filter_indel` | `false` | Filter reads with indels |

## Resource Limits

| Parameter | Default | Description |
|:----------|:--------|:------------|
| `--max_cpus` | `16` | Maximum CPUs per job |
| `--max_memory` | `128.GB` | Maximum memory per job |
| `--max_time` | `240.h` | Maximum runtime per job |

## Performance Benchmarks

Representative metrics from cfDNA duplex BAM samples:

| Sample Type | BAM Size | Variants | Runtime | CPUs |
|:------------|:---------|:---------|:--------|:-----|
| ctDNA (plasma) | 1.3 GB | 608 | ~25s | 4 |
| Plasma control | 776 MB | 608 | ~20s | 4 |

## Execution Profiles

| Profile | Description |
|:--------|:------------|
| `docker` | Local with Docker containers |
| `singularity` | HPC with Singularity |
| `slurm` | SLURM cluster with Singularity |
| `local` | No container (requires local install) |

## Related

- [Samplesheet](samplesheet.md) — Input format
- [Examples](examples.md) — Usage patterns
- [CLI Reference](../cli/run.md) — Underlying command
