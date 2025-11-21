# Advanced Usage

## Thread Configuration

Control parallelism with the `--threads` argument:

```bash
# Single-threaded (default, safest)
gbcms run ... --threads 1

# Multi-threaded for faster processing
gbcms run ... --threads 4

# Use all available cores
gbcms run ... --threads 0
```

**How it works:**
- Variants are processed in parallel
- Each thread has its own BAM reader
- Thread-safe random access to BAM files

## Custom Sample IDs

Override the default sample ID (derived from BAM filename):

```bash
# CLI argument
gbcms run ... --bam "MySampleID:/path/to/sample.bam"

# BAM list file
cat bams.txt
TumorSample    /data/tumor.bam
NormalSample   /data/normal.bam
```

Sample IDs are used in:
- Output filenames: `{sample_id}{suffix}.{format}`
- VCF sample columns
- MAF `Tumor_Sample_Barcode` field

## Output Customization

### Custom Suffixes

```bash
gbcms run ... --suffix .genotyped
# Output: sample.genotyped.vcf
```

### Output Directory Structure

```bash
gbcms run --output-dir results/batch1/
# Creates: results/batch1/{sample_id}.vcf
```

## Filter Combinations

### Strict Quality Filters

```bash
gbcms run ... \
  --min-mapping-quality 30 \
  --min-base-quality 20 \
  --filter-duplicates \
  --filter-secondary \
  --filter-supplementary \
  --filter-qc-failed
```

### Paired-End Only

```bash
gbcms run ... \
  --filter-improper-pair
```

### Exclude Reads with Indels

```bash
gbcms run ... \
  --filter-indel
```

## Docker Usage

### Basic Run

```bash
docker run -v /data:/data ghcr.io/msk-access/py-gbcms:2.0.0 \
  gbcms run \
  --fasta /data/reference.fa \
  --bam /data/sample.bam \
  --variants /data/variants.vcf \
  --output-dir /data/results/
```

### Interactive Mode

```bash
docker run -it -v /data:/data ghcr.io/msk-access/py-gbcms:2.0.0 bash
# Inside container:
gbcms run ...
```

## Performance Tips

### 1. Index Your Files

Ensure all input files are indexed:
```bash
samtools index sample.bam
samtools faidx reference.fa
tabix -p vcf variants.vcf.gz  # For compressed VCF
```

### 2. Optimal Thread Count

- **Small datasets (<10K variants)**: 1-2 threads
- **Medium datasets (10K-100K)**: 4-8 threads
- **Large datasets (>100K)**: 8-16 threads

### 3. Batch Processing

Process multiple samples in parallel using shell:

```bash
# GNU Parallel
parallel -j 4 gbcms run --fasta ref.fa --variants vars.vcf --bam {} --output-dir results/ ::: *.bam

# Simple loop
for bam in *.bam; do
  gbcms run --fasta ref.fa --variants vars.vcf --bam $bam --output-dir results/ &
done
wait
```

## Troubleshooting

### BAM Index Missing

```
Error: Could not open BAM index
```

**Solution:**
```bash
samtools index your_file.bam
```

### Chromosome Mismatch

```
Warning: BAM may not contain variant chromosomes
```

**Cause:** VCF uses `chr1` but BAM uses `1` (or vice versa)

**Solution:** Ensure consistent chromosome naming between files

### Out of Memory

**Solution:** Reduce thread count:
```bash
gbcms run ... --threads 1
```

### Slow Performance

**Checklist:**
- ✅ Files are indexed
- ✅ Using multiple threads
- ✅ BAM files are on fast storage (SSD, not NFS)
- ✅ Variant list is pre-filtered

## Integration Examples

### Nextflow Pipeline

```groovy
process countBases {
  input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path variants
  
  output:
    path "${sample_id}.vcf"
  
  script:
  """
  gbcms run \
    --fasta ${reference} \
    --bam ${bam} \
    --variants ${variants} \
    --output-dir . \
    --threads ${task.cpus}
  """
}
```

### Snakemake Workflow

```python
rule count_bases:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai",
        ref="reference.fa",
        vcf="variants.vcf"
    output:
        "results/{sample}.vcf"
    threads: 4
    shell:
        """
        gbcms run \
          --fasta {input.ref} \
          --bam {input.bam} \
          --variants {input.vcf} \
          --output-dir results/ \
          --threads {threads}
        """
```

## Environment Variables

Currently, py-gbcms does not use environment variables for configuration. All settings are passed via CLI arguments.
