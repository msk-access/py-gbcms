# Frequently Asked Questions (FAQ)

## General Questions

### What is gbcms?

gbcms is a tool for calculating base counts in BAM files at variant positions specified in VCF or MAF files. It counts how many reads support the reference allele vs. alternate allele at each position.

### Why use the Python version instead of C++?

The Python version offers:
- **Better performance**: 50-100x faster with Numba optimization
- **Type safety**: Runtime validation with Pydantic
- **Modern features**: Beautiful CLI, better error messages
- **Scalability**: Distributed computing with Ray
- **Easier to extend**: Python is more maintainable

### Is it compatible with the C++ version?

Yes! The Python version replicates 100% of C++ functionality and produces identical output.

---

## Installation

### How do I install gbcms?

```bash
# Basic installation
uv pip install gbcms

# With all features (recommended)
uv pip install "gbcms[all]"
```

### What Python version do I need?

Python 3.9 or higher. Python 3.11+ recommended for best performance.

### Do I need to install samtools?

Yes, samtools is required for creating BAM and FASTA indices:
```bash
# macOS
brew install samtools

# Ubuntu/Debian
sudo apt-get install samtools
```

### How do I install with Ray support?

```bash
uv pip install "gbcms[ray]"
```

### How do I install with fast VCF parsing?

```bash
uv pip install "gbcms[fast]"
```

---

## Usage

### What's the basic command?

```bash
gbcms count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt
```

### How do I process multiple BAM files?

```bash
# Option 1: Multiple --bam flags
gbcms count run \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --bam sample3:sample3.bam \
    ...

# Option 2: BAM file-of-files
gbcms count run --bam-fof bam_files.txt ...
```

### How do I use MAF files instead of VCF?

```bash
gbcms count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --maf variants.maf \
    --output counts.maf \
    --omaf
```

### Can I use both VCF and MAF?

No, `--vcf` and `--maf` are mutually exclusive. Choose one format.

### How do I validate files before processing?

```bash
gbcms validate files \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf
```

---

## Performance

### How can I make it faster?

1. **Install with all features**:
   ```bash
   uv pip install "gbcms[all]"
   ```

2. **Use more threads**:
   ```bash
   gbcms count run --thread 16 ...
   ```

3. **Use compressed VCF**:
   ```bash
   bgzip variants.vcf
   gbcms count run --vcf variants.vcf.gz ...
   ```

### What's the difference between counter.py and numba_counter.py?

- **counter.py**: Pure Python, flexible, baseline performance
- **numba_counter.py**: JIT-compiled, 50-100x faster, production use

Numba is used automatically when available.

### When should I use Ray?

Use Ray for:
- Very large datasets (>1M variants)
- Multi-node clusters
- Long-running jobs

Don't use Ray for:
- Small datasets (<10K variants)
- Single machine with few cores

### How much faster is cyvcf2?

cyvcf2 provides **10-100x faster** VCF loading:
- 1M variants: 195 sec → 1.8 sec (108x faster)
- Memory: 2.5 GB → 450 MB (5.5x less)

---

## Counting

### What's the difference between DMP and generic counting?

- **DMP (default)**: Specialized methods per variant type, faster
- **Generic (`--generic-counting`)**: Universal algorithm, better for complex variants

### When should I use generic counting?

Use `--generic-counting` for:
- Complex variants (MNPs, complex indels)
- Unusual ref/alt combinations
- Debugging counting discrepancies

### What are fragment counts?

Fragment counts (`--fragment-count`) count at the fragment (read pair) level instead of individual reads. Useful for detecting strand bias.

### What does fragment fractional weight do?

With `--fragment-fractional-weight`:
- Fragments with both ref and alt: RDF += 0.5, ADF += 0.5
- Without flag: Such fragments are discarded

---

## Filtering

### What quality filters are available?

```bash
--maq 20                    # Mapping quality threshold
--baq 0                     # Base quality threshold
--filter-duplicate          # Filter duplicate reads
--filter-improper-pair      # Filter improper pairs
--filter-qc-failed          # Filter QC failed reads
--filter-indel              # Filter reads with indels
--filter-non-primary        # Filter non-primary alignments
```

### What are the default filter settings?

- Mapping quality: 20
- Base quality: 0
- Filter duplicates: ON
- All other filters: OFF

### How do I disable duplicate filtering?

```bash
gbcms count run --no-filter-duplicate ...
```

---

## Output

### What output formats are supported?

1. **VCF-like** (default): Tab-delimited with counts
2. **MAF** (`--omaf`): MAF format with count columns
3. **Fillout**: Extended MAF with all samples

### What do the count columns mean?

- **DP**: Total depth
- **RD**: Reference depth
- **AD**: Alternate depth
- **DPP/RDP/ADP**: Positive strand counts
- **DPN/RDN/ADN**: Negative strand counts
- **DPF/RDF/ADF**: Fragment counts

### How do I get strand counts?

```bash
# Positive strand
gbcms count run --positive-count ...

# Negative strand
gbcms count run --negative-count ...

# Both
gbcms count run --positive-count --negative-count ...
```

---

## Troubleshooting

### Error: "command not found: gbcms"

**Solution**: Add installation directory to PATH:
```bash
export PATH="$HOME/.local/bin:$PATH"
```

### Error: "File not found"

**Check**:
- File path is correct
- File exists and is readable
- Use absolute paths if needed

### Error: "Index file not found"

**Solution**: Create indices:
```bash
samtools faidx reference.fa
samtools index sample.bam
```

### Error: "cyvcf2 not found"

**Solution**: Install with fast VCF parsing:
```bash
uv pip install "gbcms[fast]"
```

Or skip cyvcf2 (uses pure Python):
```bash
uv pip install gbcms
```

### Warning: "overlapping multimapped alignment"

**Meaning**: Same fragment end appears multiple times (multimapping)

**Action**: Usually safe to ignore. Controlled by `--suppress-warning`.

### Slow performance

**Check**:
1. Is cyvcf2 installed? (`python3 -c "import cyvcf2"`)
2. Are you using enough threads? (`--thread 16`)
3. Is VCF compressed? (`.vcf.gz` is faster)
4. Is Numba working? (Check logs for "JIT")

---

## Docker

### How do I use Docker?

```bash
# Pull image
docker pull mskaccess/gbcms:latest

# Run
docker run -v /data:/data mskaccess/gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.txt
```

### How do I build the Docker image?

```bash
docker build -t gbcms:latest .
```

---

## Advanced

### Can I use gbcms as a Python library?

Yes!

```python
from gbcms.config import Config
from gbcms.processor import VariantProcessor

config = Config(
    fasta_file="reference.fa",
    bam_files={"sample1": "sample1.bam"},
    variant_files=["variants.vcf"],
    output_file="counts.txt",
)

processor = VariantProcessor(config)
processor.process()
```

### How do I contribute?

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

### Where can I report bugs?

[GitHub Issues](https://github.com/msk-access/gbcms/issues)

### How do I get help?

1. Check this FAQ
2. Read the [documentation](README.md)
3. Search [GitHub Issues](https://github.com/msk-access/gbcms/issues)
4. Ask on [GitHub Discussions](https://github.com/msk-access/gbcms/discussions)
5. Email: access@mskcc.org

---

## Comparison with C++

### Is the Python version slower?

No! With optimizations:
- Basic: Similar speed (0.8-1.2x)
- With Numba: 2-5x faster
- With Numba + parallelization: 5-50x faster

### Does it produce identical output?

Yes, when using the same settings, output is identical.

### Can I use the same commands?

Almost! CLI flags are similar but use dashes instead of underscores:
- C++: `--filter_duplicate 1`
- Python: `--filter-duplicate`

---

## Best Practices

### What's the recommended installation?

```bash
uv pip install "gbcms[all]"
```

This includes all performance features.

### What's the recommended workflow?

1. **Validate files**:
   ```bash
   gbcms validate files ...
   ```

2. **Run with optimizations**:
   ```bash
   gbcms count run \
       --thread 16 \
       --backend joblib \
       ...
   ```

3. **Check output**:
   ```bash
   head counts.txt
   ```

### How should I organize my files?

```
project/
├── reference.fa
├── reference.fa.fai
├── bams/
│   ├── sample1.bam
│   ├── sample1.bam.bai
│   ├── sample2.bam
│   └── sample2.bam.bai
├── variants/
│   ├── variants.vcf.gz
│   └── variants.vcf.gz.tbi
└── output/
    └── counts.txt
```

---

## Quick Reference

### Installation
```bash
uv pip install "gbcms[all]"
```

### Basic Usage
```bash
gbcms count run \
    --fasta ref.fa \
    --bam s1:s1.bam \
    --vcf vars.vcf \
    --output out.txt
```

### With All Features
```bash
gbcms count run \
    --fasta ref.fa \
    --bam-fof bams.txt \
    --vcf vars.vcf.gz \
    --output out.txt \
    --thread 16 \
    --backend joblib \
    --positive-count \
    --fragment-count
```

### Validation
```bash
gbcms validate files \
    --fasta ref.fa \
    --bam s1:s1.bam \
    --vcf vars.vcf
```

---

Still have questions? Check the [full documentation](README.md) or [open an issue](https://github.com/msk-access/gbcms/issues)!
