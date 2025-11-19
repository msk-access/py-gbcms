# Input & Output Formats

## Input Formats

### Variants
gbcms supports **VCF** and **MAF** formats for defining the variants to count.

*   **VCF (`.vcf`, `.vcf.gz`)**: Standard format. Coordinates are 1-based.
*   **MAF (`.maf`)**: Mutation Annotation Format.
    *   **Requirement**: When using MAF input, you **must** provide the reference FASTA. gbcms uses the FASTA to normalize MAF indels (which often use `-`) into the VCF-style anchor-based representation required for accurate counting.

### Alignments
*   **BAM (`.bam`)**: Standard binary alignment map format. Must be indexed (`.bai`).

## Output Formats

gbcms produces one output file per input BAM sample, named `{sample_name}.{ext}`.

### VCF Output (`--format vcf`)
A standard VCF file containing the input variants with added FORMAT fields for counts.

**Header Fields:**
*   `DP`: Total Depth
*   `RD`: Reference Depth
*   `AD`: Alternate Depth
*   `VAF`: Variant Allele Frequency
*   `FSB_PVAL`: Fisher's Exact Test p-value for strand bias (fragment-level).
*   `FSB_OR`: Odds Ratio for strand bias.

### MAF Output (`--format maf`)
A tab-separated file containing standard GDC MAF columns plus custom count columns.

**Standard Columns:**
*   `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Tumor_Sample_Barcode`, `t_depth`, `t_ref_count`, `t_alt_count`.

**Custom Columns:**
*   `vcf_region`: `CHROM:POS:REF:ALT` (VCF-style representation).
*   `fsb_pval`: Fisher's Exact Test p-value for strand bias.
*   `fsb_or`: Odds Ratio for strand bias.
