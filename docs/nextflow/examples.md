# Nextflow Examples

Common usage patterns for the py-gbcms Nextflow pipeline.

## Basic Usage

### Local with Docker

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile docker
```

### SLURM Cluster

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --outdir results \
    -profile slurm
```

## MAF Output

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants input.maf \
    --fasta reference.fa \
    --format maf \
    -profile docker
```

## Strict Filtering

Enable all quality filters:

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --filter_duplicates true \
    --filter_secondary true \
    --filter_supplementary true \
    --filter_qc_failed true \
    -profile docker
```

## Resume Failed Run

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile docker \
    -resume
```

## Custom Resource Limits

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    --max_cpus 8 \
    --max_memory 32.GB \
    -profile docker
```

## Related

- [Samplesheet](samplesheet.md) — Input format
- [Parameters](parameters.md) — All options
- [CLI Quick Start](../getting-started/quickstart.md) — For few samples
