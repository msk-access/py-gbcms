# gbcms: Genomic Base Count and Analysis

**Complete orientation-aware counting system for genomic variants**

## Features

- **Complete variant analysis**: SNP, DNP, insertion, deletion, and complex variants
- **Orientation-aware counting**: Forward and reverse strand analysis
- **Fragment counting**: Proper fragment-level analysis with orientation using **Majority Rule** approach
- **Statistical analysis**: Strand bias calculation using Fisher's exact test
- **Multiple input formats**: VCF and MAF file support
- **Multiple output formats**: VCF, MAF, and Fillout format support
- **Comprehensive filtering**: 7 quality control conditions

## Installation

```bash
pip install gbcms
```

## Quick Start

```bash
# Count variants from VCF with BAM alignments
gbcms count run --fasta reference.fa --bam sample1:sample1.bam --vcf variants.vcf --output results.vcf

# Count variants from MAF with BAM alignments (sample-agnostic MAF output)
gbcms count run --fasta reference.fa --bam sample1:sample1.bam --maf variants.maf --output results.maf

# Enable fragment counting
gbcms count run --fasta reference.fa --bam sample1:sample1.bam --vcf variants.vcf --output results.vcf --fragment-count
```

## Variant Classification and Counting Strategy

gbcms automatically classifies variants and uses appropriate counting methods:

### Variant Classification
- **SNP**: Single nucleotide variants (e.g., A→T)
- **Insertion**: Insertions after single reference base (e.g., A→ATC)
- **Deletion**: Deletions to single alternate base (e.g., ATC→A)
- **Complex**: Multi-base variants and complex indels

### Counting Methods
- **Standard counting**: Optimized methods for SNP, insertion, and deletion variants
- **Generic counting**: Complete analysis for complex variants with CIGAR string parsing

All methods provide identical output including read counts, strand analysis, fragment counts, and statistical analysis readiness.

## Output Formats

- **VCF**: Standard VCF format with custom fields for counts and statistics
- **MAF**: Mutation Annotation Format with counting information
- **Fillout**: Broad Institute Firehose format compatibility
