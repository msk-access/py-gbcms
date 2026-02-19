# CLI Quick Start

Process variants with the standalone CLI.

> **Many samples on HPC?** Use [Nextflow](../nextflow/index.md) instead.

---

## Basic Usage

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/sample.vcf`

---

## Common Options

### Output Format

```bash
# VCF output (default)
gbcms run -v variants.vcf -b sample.bam -f ref.fa -o out/ --format vcf

# MAF output
gbcms run -v variants.maf -b sample.bam -f ref.fa -o out/ --format maf
```

### Multiple Samples

```bash
# Using BAM list file
echo "sample1 /path/to/sample1.bam" > bam_list.txt
echo "sample2 /path/to/sample2.bam" >> bam_list.txt

gbcms run \
    --variants variants.vcf \
    --bam-list bam_list.txt \
    --fasta reference.fa \
    --output-dir results/
```

### Custom Sample ID

```bash
gbcms run \
    --variants variants.vcf \
    --bam MySample:sample.bam \
    --fasta reference.fa \
    --output-dir results/
```
**Output:** `results/MySample.vcf`

### Quality Filters

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample.bam \
    --fasta reference.fa \
    --output-dir results/ \
    --min-mapq 30 \
    --min-baseq 20 \
    --filter-duplicates \
    --filter-secondary
```

### Threading

```bash
gbcms run ... --threads 8
```

---

## Complete Example

```bash
gbcms run \
    --variants variants.vcf \
    --bam TumorSample:tumor.bam \
    --fasta hg38.fa \
    --output-dir genotyped/ \
    --format vcf \
    --suffix .genotyped \
    --threads 8 \
    --min-mapq 30 \
    --min-baseq 20 \
    --filter-duplicates \
    --filter-secondary \
    --filter-supplementary
```

**Output:** `genotyped/TumorSample.genotyped.vcf`

---

## Docker

```bash
docker run --rm -v $(pwd):/data ghcr.io/msk-access/py-gbcms:2.7.0 \
    gbcms run \
    --variants /data/variants.vcf \
    --bam /data/sample.bam \
    --fasta /data/reference.fa \
    --output-dir /data/results/
```

---

## CLI Reference

```bash
gbcms run --help
```

| Option | Default | Description |
|:-------|:--------|:------------|
| `--variants` | Required | VCF or MAF file |
| `--bam` | Required | BAM file(s) |
| `--fasta` | Required | Reference FASTA |
| `--output-dir` | Required | Output directory |
| `--format` | vcf | Output format (vcf/maf) |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--min-baseq` | 20 | Minimum base quality |
| `--threads` | 1 | Number of threads |

---

## Next Steps

- **[Nextflow](../nextflow/index.md)** — Process many samples in parallel
- **[Architecture](../reference/architecture.md)** — How it works
- **[Allele Classification](../reference/allele-classification.md)** — How each variant type is counted
