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
| `--column-prefix` | `''` | Prefix for gbcms count columns in MAF output |
| `--preserve-barcode` | `false` | Keep original Tumor_Sample_Barcode from input MAF |
| `--show-normalization` | `false` | Append `norm_*` columns showing left-aligned coordinates |
| `--context-padding` | `5` | Flanking reference bases for haplotype construction (1-50) |
| `--threads` | `1` | Number of threads |

## Filtering Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--min-mapq` | `20` | Minimum MAPQ |
| `--min-baseq` | `20` | Minimum BASEQ |
| `--filter-duplicates` | `true` | Filter duplicate reads |
| `--filter-secondary` | `false` | Filter secondary alignments |
| `--filter-supplementary` | `false` | Filter supplementary alignments |
| `--filter-qc-failed` | `false` | Filter QC failed reads |
| `--filter-improper-pair` | `false` | Filter improperly paired reads |
| `--filter-indel` | `false` | Filter reads with indels |
| `--fragment-qual-threshold` | `10` | Quality difference threshold for fragment consensus (see [Fragment Counting](../reference/counting-metrics.md#fragment-counting)) |

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

### With Normalization Columns

```bash
gbcms run \
    --variants mutations.maf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --format maf \
    --show-normalization
```

## Output

The command produces a VCF or MAF file with:

- **Allele counts** (reference and alternate depth)
- **VAF** (variant allele frequency)
- **Strand bias** (Fisher's exact test)
- **Fragment counts** (deduplicated)
- **Validation status** (`PASS`, `PASS_WARN_HOMOPOLYMER_DECOMP`, `REF_MISMATCH`, or `FETCH_FAILED`)
- **Normalization columns** (with `--show-normalization`): left-aligned position, REF, and ALT

## Related

- [Quick Start](../getting-started/quickstart.md) — Common patterns
- [gbcms normalize](normalize.md) — Standalone normalization (no counting)
- [Nextflow Pipeline](../nextflow/index.md) — For many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specs
- [Variant Counting](../reference/allele-classification.md) — How each variant type is counted
