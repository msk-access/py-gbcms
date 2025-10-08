# Input & Output Formats

Complete guide to input and output formats supported by gbcms.

## Input Formats

### VCF Files

**Format**: Variant Call Format (VCF)

**Supported versions**: VCF 4.0, 4.1, 4.2

**File extensions**:
- `.vcf` - Uncompressed VCF
- `.vcf.gz` - Compressed VCF (with cyvcf2)

**Required columns**:
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    100  .   A    T    .     PASS    .
```

**Usage**:
```bash
gbcms count run --vcf variants.vcf ...
gbcms count run --vcf variants.vcf.gz ...  # Faster with cyvcf2
```

**Multiple VCF files**:
```bash
gbcms count run \
    --vcf variants1.vcf \
    --vcf variants2.vcf \
    --vcf variants3.vcf \
    ...
```

### MAF Files

**Format**: Mutation Annotation Format (MAF)

**File extension**: `.maf`

**Required columns**:
- `Hugo_Symbol`
- `Chromosome`
- `Start_Position`
- `End_Position`
- `Reference_Allele`
- `Tumor_Seq_Allele1`
- `Tumor_Seq_Allele2`
- `Tumor_Sample_Barcode`
- `Matched_Norm_Sample_Barcode`
- `Variant_Classification`

**Usage**:
```bash
gbcms count run --maf variants.maf ...
```

**Note**: `--maf` and `--vcf` are mutually exclusive.

### BAM Files

**Format**: Binary Alignment Map (BAM)

**Requirements**:
- Must be coordinate-sorted
- Must have index file (`.bai`)

**Single BAM**:
```bash
gbcms count run --bam sample1:sample1.bam ...
```

**Multiple BAMs**:
```bash
gbcms count run \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --bam sample3:sample3.bam \
    ...
```

**BAM File-of-Files**:
```bash
# Create bam_files.txt
cat > bam_files.txt << EOF
sample1	sample1.bam
sample2	sample2.bam
sample3	sample3.bam
EOF

# Use file-of-files
gbcms count run --bam-fof bam_files.txt ...
```

### Reference FASTA

**Format**: FASTA format

**Requirements**:
- Must have index file (`.fai`)

**Usage**:
```bash
gbcms count run --fasta reference.fa ...
```

**Create index**:
```bash
samtools faidx reference.fa
```

---

## Output Formats

### VCF Format (Proper VCF with INFO fields)

**Extension**: `.vcf`

**Structure**: Standard VCF format with count and strand bias information in FORMAT and INFO fields

**Example**:
```vcf
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">
##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth for this sample">
##FORMAT=<ID=SB,Number=3,Type=Float,Description="Strand bias for this sample">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO                    FORMAT      SAMPLE1             SAMPLE2
chr1    100  .   A    T    .     .       DP=95;SB=0.001234,2.5,reverse  DP:RD:AD:SB  50:30:20:0.001234:2.5:reverse  45:25:20:0.95:1.2:none
```

**INFO Fields**:
- `DP`: Total depth across all samples
- `SB`: Strand bias (p-value,odds_ratio,direction) - most significant across samples
- `FSB`: Fragment strand bias (when fragment counting enabled)

**FORMAT Fields**:
- `DP`: Total depth for this sample
- `RD`: Reference allele depth for this sample
- `AD`: Alternate allele depth for this sample
- `DPP`: Positive strand depth (if enabled)
- `RDP`: Positive strand reference depth (if enabled)
- `ADP`: Positive strand alternate depth (if enabled)
- `DPF`: Fragment depth (if enabled)
- `RDF`: Fragment reference depth (if enabled)
- `ADF`: Fragment alternate depth (if enabled)
- `SB`: Strand bias (p-value,odds_ratio,direction) for this sample
- `FSB`: Fragment strand bias for this sample (if enabled)

**Usage**:
```bash
gbcms count run --output counts.vcf ...
```

**Compatible Tools**:
- ✅ `bcftools view`, `bcftools filter`
- ✅ GATK tools
- ✅ VCF parsing libraries

### MAF vs Fillout Format Comparison

| **Aspect** | **MAF Format** | **Fillout Format** |
|------------|----------------|-------------------|
| **Best For** | Tumor-normal studies | Multi-sample cohorts |
| **Structure** | Tumor-Normal pairs | All samples as columns |
| **Layout** | Long format (t_*, n_*) | Wide format (sample:DP) |
| **Compatibility** | TCGA standard | gbcms-specific |
| **Use Case** | \`--omaf\` flag | Default with MAF input |

#### When to Use MAF Format
- ✅ **Tumor-Normal Studies**: Matched tumor-normal sample pairs
- ✅ **TCGA Integration**: Compatible with TCGA tools and standards
- ✅ **Standard Workflows**: Traditional cancer genomics analysis

#### When to Use Fillout Format
- ✅ **Cohort Studies**: Multiple samples analyzed together
- ✅ **Population Genomics**: Large-scale sample collections
- ✅ **Research Pipelines**: Complex multi-sample analysis workflows

### MAF Format (Tumor-Normal Structure)

**Extension**: `.maf`

**Structure**: Standard TCGA-compatible MAF with tumor-normal columns

**Example Output**:
```tsv
Hugo_Symbol  Chromosome  Start_Position  t_depth  t_ref_count  t_alt_count  n_depth  n_ref_count  n_alt_count
TP53        chr17       7577120        150      75           75           120      110          10
```

**Usage**:
```bash
gbcms count run --fasta ref.fa --bam tumor.bam --maf variants.maf --omaf --output results.maf
```

### Fillout Format (Multi-Sample Structure)

**Extension**: `.txt` or `.maf`

**Structure**: Extended MAF with all samples represented as columns

**Example Output**:
```tsv
Hugo_Symbol  sample1_DP  sample1_RD  sample1_AD  sample2_DP  sample2_RD  sample2_AD
TP53        150         75          75          200         180         20
```

**Usage**:
```bash
gbcms count run --fasta ref.fa --bam sample_*.bam --maf variants.maf --output results.fillout
```

---

## Count Types

### Read-Level Counts

**DP** (Depth):
- Total number of reads covering the position
- Includes both reference and alternate reads

**RD** (Reference Depth):
- Number of reads matching the reference allele

**AD** (Alternate Depth):
- Number of reads matching the alternate allele

### Strand Counts

**Positive Strand** (enabled with `--positive-count`):
- `DPP`: Positive strand depth
- `RDP`: Positive strand reference depth
- `ADP`: Positive strand alternate depth

**Negative Strand** (enabled with `--negative-count`):
- `DPN`: Negative strand depth
- `RDN`: Negative strand reference depth
- `ADN`: Negative strand alternate depth

### Fragment Counts

**Fragment-Level** (enabled with `--fragment-count`):
- `DPF`: Fragment depth (number of unique fragments)
- `RDF`: Fragments with reference allele
- `ADF`: Fragments with alternate allele

**Fractional Weights** (with `--fragment-fractional-weight`):
- Fragments with both ref and alt: RDF += 0.5, ADF += 0.5
- Fragments with only ref: RDF += 1
- Fragments with only alt: ADF += 1

### Strand Bias Analysis ⭐

**Statistical Strand Bias** (automatically calculated):
- `SB_PVAL`: Strand bias p-value from Fisher's exact test
- `SB_OR`: Strand bias odds ratio
- `SB_DIR`: Strand bias direction ("forward", "reverse", "none")
- `FSB_PVAL`: Fragment strand bias p-value (if fragment counting enabled)
- `FSB_OR`: Fragment strand bias odds ratio (if fragment counting enabled)
- `FSB_DIR`: Fragment strand bias direction (if fragment counting enabled)

**Methodology**:
- Uses **Fisher's exact test** for statistical rigor
- **10% threshold** for direction determination
- **Minimum 10 reads** for reliable calculations
- **Automatic calculation** for all variants and samples

---

## Examples

### Example 1: VCF to VCF Output

```bash
gbcms count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --vcf variants.vcf \
    --output counts.vcf
```

**Output** (`counts.vcf`):
```vcf
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth across all samples">
##INFO=<ID=SB,Number=3,Type=Float,Description="Strand bias p-value, odds ratio, direction">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth for this sample">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference allele depth for this sample">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternate allele depth for this sample">
##FORMAT=<ID=SB,Number=3,Type=Float,Description="Strand bias for this sample">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO                    FORMAT        SAMPLE1               SAMPLE2
chr1    100  .   A    T    .     .       DP=95;SB=0.001234,2.5,reverse  DP:RD:AD:SB  50:30:20:0.001234:2.5:reverse  45:25:20:0.95:1.2:none
chr1    200  .   C    G    .     .       DP=115;SB=0.95,1.2,none       DP:RD:AD:SB  60:40:20:0.95:1.2:none        55:35:20:0.95:1.2:none
```

### Example 2: MAF to MAF Output

```bash
gbcms count run \
    --fasta reference.fa \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --maf variants.maf \
    --output counts.maf \
    --omaf
```

**Output** (`counts.maf`): Original MAF columns + count columns + strand bias columns

### Example 3: For Tabular Output, Use MAF Format

```bash
# For tabular output, use MAF format (which is naturally tabular)
gbcms count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --maf variants.maf \
    --output counts.txt \
    --omaf
```

**Note**: For tabular output, use MAF format with `--omaf`. VCF format is now proper VCF with INFO/FORMAT fields.

### Example 4: Fillout for Multiple Samples

```bash
# Create BAM file-of-files with all samples
cat > all_samples.txt << EOF
sample1	sample1.bam
sample2	sample2.bam
sample3	sample3.bam
sample4	sample4.bam
EOF

# Run fillout
gbcms count run \
    --fasta reference.fa \
    --bam-fof all_samples.txt \
    --maf somatic_variants.maf \
    --output fillout.maf \
    --omaf
```

**Result**: Counts for ALL samples at each variant position

---

## File Validation

### Validate Before Processing

```bash
gbcms validate files \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf
```

**Checks**:
- ✅ Files exist
- ✅ FASTA index exists
- ✅ BAM index exists
- ✅ Files are readable
- ✅ Formats are valid

---

## Best Practices

### 1. Use Compressed VCF

```bash
# Compress VCF
bgzip variants.vcf
tabix -p vcf variants.vcf.gz

# Use compressed (faster with cyvcf2)
gbcms count run --vcf variants.vcf.gz ...
```

### 2. Index All Files

```bash
# Index FASTA
samtools faidx reference.fa

# Index BAM
samtools index sample.bam

# Index VCF
tabix -p vcf variants.vcf.gz
```

### 3. Use BAM File-of-Files for Many Samples

```bash
# Instead of many --bam flags
gbcms count run --bam-fof samples.txt ...
```

### 4. Validate First

```bash
# Check files before long processing
gbcms validate files ...
```

---

## Troubleshooting

### Issue: "File not found"

**Check**:
- File path is correct
- File exists
- File is readable

### Issue: "Index file not found"

**Solution**:
```bash
# Create FASTA index
samtools faidx reference.fa

# Create BAM index
samtools index sample.bam
```

### Issue: "Invalid VCF format"

**Solution**:
- Check VCF has required columns
- Validate with: `bcftools view -h variants.vcf`
- Use cyvcf2 for better error handling

### Issue: "MAF missing required columns"

**Check**: Ensure MAF has all required columns listed above

---

## Summary

### Input Formats
- ✅ VCF (`.vcf`, `.vcf.gz`)
- ✅ MAF (`.maf`)
- ✅ BAM (`.bam` + `.bai`)
- ✅ FASTA (`.fa` + `.fai`)

### Output Formats
- ✅ **VCF** (proper VCF format with INFO/FORMAT fields)
- ✅ MAF (`--omaf`)
- ✅ Fillout (extended MAF)

### Count Types
- ✅ Read-level (DP, RD, AD)
- ✅ Strand counts (DPP, RDP, ADP, DPN, RDN, ADN)
- ✅ Fragment counts (DPF, RDF, ADF)
- ✅ **Strand bias analysis** (SB_PVAL, SB_OR, SB_DIR, FSB_*) ⭐

See [CLI Features](CLI_FEATURES.md) for complete command reference.
