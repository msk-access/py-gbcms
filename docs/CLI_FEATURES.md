# CLI Reference

The primary entry point for gbcms is the `run` command.

```bash
python -m gbcms.cli run [OPTIONS]
```

## Arguments

| Option | Short | Description | Required |
| :--- | :--- | :--- | :--- |
| `--variants` | `-v` | Path to VCF or MAF file containing variants to count. | **Yes** |
| `--fasta` | `-f` | Path to the reference genome FASTA file. | **Yes** |
| `--output-dir` | `-o` | Directory where output files will be written. | **Yes** |
| `--bam` | `-b` | Path to a BAM file. Can be used multiple times. | No* |
| `--bam-list` | `-L` | Path to a file containing a list of BAM paths (one per line). | No* |
| `--format` | | Output format: `vcf` or `maf`. Default: `vcf`. | No |
| `--min-mapq` | | Minimum mapping quality (MAPQ) to include a read. Default: `20`. | No |
| `--min-baseq` | | Minimum base quality to count a base. Default: `0`. | No |

*\* At least one of `--bam` or `--bam-list` must be provided.*

## Examples

### 1. Single Sample, VCF Output

```bash
python -m gbcms.cli run \
    --variants input.vcf \
    --fasta reference.fa \
    --bam sample1.bam \
    --output-dir results \
    --format vcf
```

### 2. Multiple Samples, MAF Output

```bash
python -m gbcms.cli run \
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
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bam
```

Then run:

```bash
python -m gbcms.cli run \
    --variants input.vcf \
    --fasta reference.fa \
    --bam-list bams.txt \
    --output-dir results
```
