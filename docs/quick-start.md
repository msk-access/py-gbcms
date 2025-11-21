# Quick Start

## Basic Usage

Get started with py-gbcms in minutes.

### 1. Installation

```bash
pip install py-gbcms
```

### 2. Prepare Your Data

You'll need:
- **Reference FASTA**: Indexed genome reference (`.fa` + `.fa.fai`)
- **BAM File(s)**: Aligned reads (`.bam` + `.bam.bai`)
- **Variants**: VCF or MAF file with variants to count

### 3. Run Your First Command

**Single BAM, VCF output:**

```bash
gbcms run \
  --fasta reference.fa \
  --bam sample.bam \
  --variants variants.vcf \
  --output-dir results/
```

**Multiple BAMs, MAF output:**

```bash
gbcms run \
  --fasta reference.fa \
  --bam tumor.bam \
  --bam normal.bam \
  --variants variants.maf \
  --format maf \
  --output-dir results/
```

### 4. Check Your Results

Output files are named: `{sample_id}{suffix}.{format}`

```bash
ls results/
# sample.vcf
# tumor.maf
# normal.maf
```

## Common Patterns

### Processing Multiple Samples

Create a BAM list file (`bams.txt`):

```text
SampleA  /path/to/sampleA.bam
SampleB  /path/to/sampleB.bam
SampleC  /path/to/sampleC.bam
```

Run:

```bash
gbcms run \
  --fasta reference.fa \
  --bam-list bams.txt \
  --variants variants.vcf \
  --output-dir results/
```

### Filtering Options

```bash
gbcms run ... \
  --min-mapping-quality 20 \
  --min-base-quality 10 \
  --filter-duplicates \
  --filter-secondary \
  --filter-supplementary
```

### Performance Tuning

```bash
# Use 4 threads for faster processing
gbcms run ... --threads 4

# Custom output naming
gbcms run ... --suffix .counted
```

## Next Steps

- Learn about all [CLI options](CLI_FEATURES.md)
- Understand [input/output formats](INPUT_OUTPUT.md)
- Explore [advanced usage](advanced-usage.md)
