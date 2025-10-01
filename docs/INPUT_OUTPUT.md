# Input & Output Formats

Complete guide to input and output formats supported by GetBaseCounts.

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
getbasecounts count run --vcf variants.vcf ...
getbasecounts count run --vcf variants.vcf.gz ...  # Faster with cyvcf2
```

**Multiple VCF files**:
```bash
getbasecounts count run \
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
getbasecounts count run --maf variants.maf ...
```

**Note**: `--maf` and `--vcf` are mutually exclusive.

### BAM Files

**Format**: Binary Alignment Map (BAM)

**Requirements**:
- Must be coordinate-sorted
- Must have index file (`.bai`)

**Single BAM**:
```bash
getbasecounts count run --bam sample1:sample1.bam ...
```

**Multiple BAMs**:
```bash
getbasecounts count run \
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
getbasecounts count run --bam-fof bam_files.txt ...
```

### Reference FASTA

**Format**: FASTA format

**Requirements**:
- Must have index file (`.fai`)

**Usage**:
```bash
getbasecounts count run --fasta reference.fa ...
```

**Create index**:
```bash
samtools faidx reference.fa
```

---

## Output Formats

### VCF-like Format (Default)

**Extension**: `.txt`

**Columns**:
```
Chrom  Pos  Ref  Alt  Sample1_DP  Sample1_RD  Sample1_AD  Sample1_DPP  Sample1_RDP  Sample1_ADP  Sample2_DP  ...
chr1   100  A    T    50          30          20          25           15           10           45          ...
```

**Column descriptions**:
- `DP`: Total depth
- `RD`: Reference depth
- `AD`: Alternate depth
- `DPP`: Positive strand depth
- `RDP`: Positive strand reference depth
- `ADP`: Positive strand alternate depth

**Usage**:
```bash
getbasecounts count run --output counts.txt ...
```

### MAF Format

**Extension**: `.maf`

**Columns**: All original MAF columns plus count columns

**Usage**:
```bash
getbasecounts count run --maf variants.maf --output counts.maf --omaf
```

**Count columns added**:
- `{sample}_DP`, `{sample}_RD`, `{sample}_AD`
- `{sample}_DPP`, `{sample}_RDP`, `{sample}_ADP` (if `--positive-count`)
- `{sample}_DPN`, `{sample}_RDN`, `{sample}_ADN` (if `--negative-count`)
- `{sample}_DPF`, `{sample}_RDF`, `{sample}_ADF` (if `--fragment-count`)

### Fillout Format

**Extension**: `.txt` or `.maf`

**Description**: Extended MAF format with counts for ALL samples

**Usage**:
```bash
getbasecounts count run \
    --maf variants.maf \
    --bam-fof all_samples.txt \
    --output fillout.txt \
    --omaf
```

**Use case**: Genotyping all samples at variant positions

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

---

## Examples

### Example 1: VCF to VCF-like Output

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --vcf variants.vcf \
    --output counts.txt
```

**Output** (`counts.txt`):
```
Chrom  Pos  Ref  Alt  sample1_DP  sample1_RD  sample1_AD  sample2_DP  sample2_RD  sample2_AD
chr1   100  A    T    50          30          20          45          25          20
chr1   200  C    G    60          40          20          55          35          20
```

### Example 2: MAF to MAF Output

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --maf variants.maf \
    --output counts.maf \
    --omaf
```

**Output** (`counts.maf`): Original MAF columns + count columns

### Example 3: With All Count Types

```bash
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts_full.txt \
    --positive-count \
    --negative-count \
    --fragment-count
```

**Output columns**:
```
Chrom  Pos  Ref  Alt  
sample1_DP  sample1_RD  sample1_AD  
sample1_DPP sample1_RDP sample1_ADP  
sample1_DPN sample1_RDN sample1_ADN  
sample1_DPF sample1_RDF sample1_ADF
```

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
getbasecounts count run \
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
getbasecounts validate files \
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
getbasecounts count run --vcf variants.vcf.gz ...
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
getbasecounts count run --bam-fof samples.txt ...
```

### 4. Validate First

```bash
# Check files before long processing
getbasecounts validate files ...
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
- ✅ VCF-like (default)
- ✅ MAF (`--omaf`)
- ✅ Fillout (extended MAF)

### Count Types
- ✅ Read-level (DP, RD, AD)
- ✅ Strand counts (DPP, RDP, ADP, DPN, RDN, ADN)
- ✅ Fragment counts (DPF, RDF, ADF)

See [CLI Features](CLI_FEATURES.md) for complete command reference.
