# gbcms run

Count alleles at variant positions across one or more BAM files.

## Synopsis

```bash
gbcms run [OPTIONS] --variants <FILE> --bam <NAME:PATH>... --fasta <FILE>
```

## Required Arguments

| Option | Description |
|:-------|:------------|
| `--variants` | [VCF or MAF](../reference/input-formats.md) file with variant positions |
| `--bam` | BAM file in `name:path` format (can repeat) |
| `--fasta` | Reference FASTA file (with .fai index) |

## Output Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--output-dir` | `.` | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |
| `--threads` | `4` | Number of threads |

## Filtering Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--min-mapq` | `20` | Minimum MAPQ |
| `--min-baseq` | `0` | Minimum BASEQ |
| `--filter-duplicates` | `false` | Filter duplicate reads |
| `--filter-secondary` | `false` | Filter secondary alignments |
| `--filter-supplementary` | `false` | Filter supplementary alignments |
| `--filter-qc-failed` | `false` | Filter QC failed reads |
| `--filter-improper-pair` | `false` | Filter improperly paired reads |
| `--filter-indel` | `false` | Filter reads with indels |

## Examples

### Single BAM

```bash
gbcms run \
    --variants mutations.vcf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --output-dir results/
```

### Multiple BAMs

```bash
gbcms run \
    --variants mutations.maf \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --fasta reference.fa \
    --format maf
```

### With Filtering

```bash
gbcms run \
    --variants mutations.vcf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --filter-duplicates \
    --min-mapq 30
```

## Output

The command produces a VCF or MAF file with:

- **Allele counts** (reference and alternate depth)
- **VAF** (variant allele frequency)
- **Strand bias** (Fisher's exact test)
- **Fragment counts** (deduplicated)

## Related

- [Quick Start](../getting-started/quickstart.md) — Common patterns
- [Nextflow Pipeline](../nextflow/index.md) — For many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specs
