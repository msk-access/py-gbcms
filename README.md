# gbcms: Genomic Base Count and Analysis

**Complete orientation-aware counting system for genomic variants**

## Features

- **High Performance**: **Rust-powered** core engine with multi-threading support for maximum speed.
- **Complete variant analysis**: SNP, MNP, insertion, deletion, and **complex variants** (DelIns, SNP+Indel)
- **Orientation-aware counting**: Forward and reverse strand analysis
- **Fragment counting**: Proper fragment-level analysis with orientation using **Majority Rule** approach
- **Statistical analysis**: Strand bias calculation using Fisher's exact test
- **Multiple input formats**: VCF and MAF file support
- **Multiple output formats**: VCF and MAF support
- **Comprehensive filtering**: 7 quality control conditions

## Installation

```bash
pip install gbcms
```

## Quick Start

```bash
# Count variants from VCF with BAM alignments
python -m gbcms.cli run --fasta reference.fa --bam sample1.bam --variants variants.vcf --output-dir results/

# Count variants from MAF with BAM alignments (sample-agnostic MAF output)
python -m gbcms.cli run --fasta reference.fa --bam sample1.bam --variants variants.maf --output-dir results/ --format maf

# Process multiple BAMs using a File of Files (FoF)
python -m gbcms.cli run --fasta reference.fa --bam-list bam_list.txt --variants variants.vcf --output-dir results/

# Enable fragment counting (implicit in v2, always calculated)
```

## Variant Classification and Counting Strategy

gbcms automatically classifies variants and uses appropriate counting methods:

### Variant Classification
- **SNP**: Single nucleotide variants (e.g., A→T)
- **MNP**: Multi-nucleotide polymorphisms (e.g., AT→CG)
- **Insertion**: Insertions after single reference base (e.g., A→ATC)
- **Deletion**: Deletions to single alternate base (e.g., ATC→A)
- **Complex**: Multi-base variants and complex indels

### Counting Methods
- **Standard counting**: Optimized methods for SNP, insertion, and deletion variants
- **Haplotype Reconstruction**: Advanced sequence reconstruction for complex variants (e.g. DelIns, SNP+Indel) to ensure accurate counting regardless of CIGAR complexity.

All methods provide identical output including read counts, strand analysis, fragment counts, and statistical analysis readiness.

## Output Formats

- **VCF**: Standard VCF format with custom fields for counts and statistics
- **MAF**: Mutation Annotation Format with counting information
- **Fillout**: Broad Institute Firehose format compatibility
