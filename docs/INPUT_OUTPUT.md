# Input & Output Formats

## Input Formats

### Variants
gbcms supports **VCF** and **MAF** formats for defining the variants to count. It handles:
*   **SNPs**: Single Nucleotide Polymorphisms.
*   **MNPs**: Multi-Nucleotide Polymorphisms.
*   **Insertions**: Pure insertions.
*   **Deletions**: Pure deletions.
*   **Complex**: Mixed variants (e.g., DelIns, SNP+Indel). These are handled via **Haplotype Reconstruction**, where the read's sequence is reconstructed and compared to the ALT allele.

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
*   `DP`: Total Depth (`ref_total,alt_total`)
*   `RD`: Reference Read Depth (`fwd,rev`)
*   `AD`: Alternate Read Depth (`fwd,rev`)
*   `RDF`: Reference Fragment Depth (`fwd,rev`)
*   `ADF`: Alternate Fragment Depth (`fwd,rev`)
*   `VAF`: Variant Allele Fraction (Read Level)
*   `FAF`: Variant Allele Fraction (Fragment Level)
*   `SB_PVAL`: Fisher's Exact Test p-value for strand bias (read-level).
*   `SB_OR`: Odds Ratio for strand bias (read-level).
*   `FSB_PVAL`: Fisher's Exact Test p-value for strand bias (fragment-level).
*   `FSB_OR`: Odds Ratio for strand bias (fragment-level).

### MAF Output (`--format maf`)
A tab-separated file containing standard GDC MAF columns plus custom count columns.

**Standard Columns:**
*   `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Tumor_Sample_Barcode`, `t_depth`, `t_ref_count`, `t_alt_count`.

**Custom Columns:**
*   `t_ref_count_forward`, `t_ref_count_reverse`
*   `t_alt_count_forward`, `t_alt_count_reverse`
*   `t_ref_count_fragment_forward`, `t_ref_count_fragment_reverse`
*   `t_alt_count_fragment_forward`, `t_alt_count_fragment_reverse`
*   `t_vaf`: Variant Allele Fraction (Read Level)
*   `t_vaf_fragment`: Variant Allele Fraction (Fragment Level)
*   `vcf_region`: `CHROM:POS:REF:ALT` (VCF-style representation).
*   `fsb_pval`: Fisher's Exact Test p-value for strand bias.
*   `fsb_or`: Odds Ratio for strand bias.
