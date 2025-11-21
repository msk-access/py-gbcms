# CLI Reference

The primary entry point for gbcms is the `run` command.

```bash
gbcms run [OPTIONS]
```

## Arguments

| Option | Short | Description | Required |
| :--- | :--- | :--- | :--- |
| `--variants` | `-v` | Path to VCF or MAF file containing variants to count. | **Yes** |
| `--fasta` | `-f` | Path to the reference genome FASTA file. | **Yes** |
| `--output-dir` | `-o` | Directory where output files will be written. | **Yes** |
| `--bam` | `-b` | Path to a BAM file. Can be used multiple times. Supports `ID:path` format. | No* |
| `--bam-list` | `-L` | Path to a file containing a list of BAM paths (one per line, optional ID column). | No* |
| `--format` | | Output format: `vcf` or `maf`. Default: `vcf`. | No |
| `--suffix` | `-S` | Suffix to append to output filename (e.g. `.genotyped`). | No |
| `--threads` | `-t` | Number of threads. Default: `1`. | No |
| `--min-mapq` | | Minimum mapping quality (MAPQ) to include a read. Default: `20`. | No |
| `--min-baseq` | | Minimum base quality to count a base. Default: `0`. | No |
| `--filter-duplicates` | | Filter duplicate reads. Default: `True`. | No |
| `--filter-qc-failed` | | Filter reads failing platform/vendor quality checks. Default: `False`. | No |
| `--filter-improper-pair` | | Filter reads that are not marked as proper pairs. Default: `False`. | No |
| `--filter-indel` | | Filter reads containing insertions or deletions. Default: `False`. | No |
| `--filter-secondary` | | Filter secondary alignments. Default: `False`. | No |
| `--filter-supplementary` | | Filter supplementary alignments. Default: `False`. | No |

*\* At least one of `--bam` or `--bam-list` must be provided.*

## Examples

### 1. Single Sample, VCF Output

```bash
gbcms run \
    --variants input.vcf \
    --fasta reference.fa \
    --bam sample1.bam \
    --output-dir results \
    --format vcf
```

### 2. Multiple Samples, MAF Output

```bash
gbcms run \
    --variants input.maf \
    --fasta reference.fa \
    --bam sample1.bam \
    --bam sample2.bam \
    --output-dir results \
    --format maf
```

### 3. File of Files (FoF)

If you have many BAM files, list them in a text file (e.g., `bams.txt`):

```text
SampleID1   /path/to/sample1.bam
SampleID2   /path/to/sample2.bam
```

Then run:

```bash
gbcms run \
    --variants input.vcf \
    --fasta reference.fa \
    --bam-list bams.txt \
    --output-dir results
```
