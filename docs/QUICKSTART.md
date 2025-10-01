# GetBaseCounts - Quick Start Guide

Get up and running with GetBaseCounts in 5 minutes!

## Installation

### Option 1: Using uv (Recommended)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install getbasecounts
uv pip install getbasecounts
```

### Option 2: Using pip

```bash
pip install getbasecounts
```

### Option 3: Using Docker

```bash
docker pull mskaccess/getbasecounts:latest
```

## Verify Installation

```bash
getbasecounts version
```

You should see a nice panel with version information.

## Basic Usage

### 1. Prepare Your Data

You'll need:
- **Reference FASTA** (indexed with `.fai`)
- **BAM files** (indexed with `.bai`)
- **Variant file** (VCF or MAF format)

### 2. Index Your Files (if not already done)

```bash
# Index reference FASTA
samtools faidx reference.fa

# Index BAM files
samtools index sample1.bam
samtools index sample2.bam
```

### 3. Validate Your Files (Optional but Recommended)

```bash
getbasecounts validate files \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf
```

This will check if all files exist and have required indices.

### 4. Run GetBaseCounts

#### Simple VCF Example

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt
```

#### Multiple Samples

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --bam sample3:sample3.bam \
    --vcf variants.vcf \
    --output counts.txt
```

#### Using BAM File-of-Files

Create `bam_files.txt`:
```
sample1	/path/to/sample1.bam
sample2	/path/to/sample2.bam
sample3	/path/to/sample3.bam
```

Then run:
```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --vcf variants.vcf \
    --output counts.txt
```

#### MAF Input/Output

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --maf variants.maf \
    --output counts.maf \
    --omaf
```

## Common Options

### Quality Filters

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --maq 30 \                      # Mapping quality threshold
    --baq 20 \                      # Base quality threshold
    --filter-duplicate \            # Filter duplicate reads
    --filter-improper-pair          # Filter improper pairs
```

### Performance Options

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 8                      # Use 8 threads
```

### Output Options

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --positive-count \              # Include positive strand counts
    --fragment-count                # Include fragment counts
```

## Output Format

### VCF Output

Tab-separated file with columns:
```
Chrom  Pos  Ref  Alt  sample1:DP  sample1:RD  sample1:AD  ...
chr1   100  A    T    50          30          20          ...
```

Where:
- **DP**: Total depth
- **RD**: Reference depth
- **AD**: Alternate depth
- **DPP**: Positive strand depth (if `--positive-count`)
- **RDP**: Positive strand reference depth
- **ADP**: Positive strand alternate depth
- **DPF**: Fragment depth (if `--fragment-count`)
- **RDF**: Fragment reference depth
- **ADF**: Fragment alternate depth

### MAF Output

Standard MAF format with updated count columns:
- `t_depth`, `t_ref_count`, `t_alt_count`
- `n_depth`, `n_ref_count`, `n_alt_count`

## Docker Usage

### Basic Run

```bash
docker run -v /path/to/data:/data getbasecounts:latest \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/output.txt
```

### Using Docker Compose

Create `docker-compose.yml`:
```yaml
version: '3.8'
services:
  getbasecounts:
    image: getbasecounts:latest
    volumes:
      - ./data:/data
    command: >
      --fasta /data/reference.fa
      --bam sample1:/data/sample1.bam
      --vcf /data/variants.vcf
      --output /data/output.txt
```

Run:
```bash
docker-compose up
```

## Troubleshooting

### Error: "Reference FASTA index not found"

**Solution**: Index your FASTA file
```bash
samtools faidx reference.fa
```

### Error: "BAM index not found"

**Solution**: Index your BAM files
```bash
samtools index sample.bam
```

### Error: "Could not find variant chrom in BAM file"

**Solution**: Ensure chromosome names match between VCF/MAF and BAM files
- VCF uses `chr1` but BAM uses `1`, or vice versa
- Check with: `samtools view -H sample.bam | grep @SQ`

### Slow Performance

**Solutions**:
1. Increase threads: `--thread 8`
2. Reduce block size for memory-constrained systems: `--max-block-size 5000`
3. Use Docker for consistent performance

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [CONTRIBUTING.md](CONTRIBUTING.md) to contribute
- Report issues on [GitHub](https://github.com/msk-access/getbasecounts/issues)

## Example Workflow

Complete example from start to finish:

```bash
# 1. Install
uv pip install getbasecounts

# 2. Prepare data
samtools faidx reference.fa
samtools index tumor.bam
samtools index normal.bam

# 3. Create BAM file list
cat > bam_files.txt << EOF
tumor	tumor.bam
normal	normal.bam
EOF

# 4. Validate files (optional)
getbasecounts validate files \
    --fasta reference.fa \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --vcf variants.vcf

# 5. Run analysis
getbasecounts count run \
    --fasta reference.fa \
    --bam-fof bam_files.txt \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 4 \
    --maq 20 \
    --baq 10 \
    --filter-duplicate

# 5. View results
head counts.txt
```

## Getting Help

```bash
# Show all commands
getbasecounts --help

# Show help for count command
getbasecounts count --help

# Show help for count run
getbasecounts count run --help

# Show help for validate
getbasecounts validate --help

# Show version
getbasecounts version

# Show tool info
getbasecounts info
```

## CLI Features

GetBaseCounts uses **Typer** with **Rich** for a beautiful CLI experience:

- 📋 **Organized help panels** - Options grouped by category
- 🎨 **Colored output** - Easy to read terminal output
- 📊 **Progress bars** - Visual feedback during processing
- ✅ **Validation** - File validation before processing
- 🔍 **Subcommands** - Clean command structure

Happy counting! 🧬
