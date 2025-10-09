# gbcms Architecture

## Overview

gbcms provides complete orientation-aware counting for genomic variants with a focus on consistency and accuracy across all variant types.

## Core Components

### Variant Classification
- **Length-based classification**: Uses reference and alternate allele lengths
- **Categories**: SNP, Insertion, Deletion, Complex
- **Format agnostic**: Same classification for VCF and MAF inputs

### Counting Strategy
- **Standard counting**: Optimized methods for simple variant types
  - `count_bases_snp()`: For SNP variants
  - `count_bases_indel()`: For insertion and deletion variants
- **Generic counting**: Complete analysis for complex variants
  - `count_bases_generic()`: For multi-base and complex variants
  - Includes CIGAR string parsing for complex indel handling

### Key Features
- **Complete strand analysis**: Forward and reverse strand counting for all alleles
- **Fragment counting**: Proper fragment-level analysis with orientation tracking
- **Statistical analysis**: Strand bias calculation using Fisher's exact test
- **Comprehensive filtering**: 7 quality control conditions applied consistently

## Data Flow

1. **Input parsing**: VCF or MAF files parsed with appropriate format handling
2. **Variant classification**: Automatic classification based on allele lengths
3. **Alignment processing**: BAM file alignment analysis with filtering
4. **Counting**: Route to appropriate counting method based on variant type
5. **Statistical analysis**: Strand bias calculation for all variant types
6. **Output formatting**: Generate VCF, MAF, or Fillout format outputs

## Consistency Guarantees

- **Identical analysis**: All variant types receive same level of analysis
- **Complete feature parity**: All counting and analysis features available for all types
- **Format consistency**: Same behavior for VCF and MAF inputs
- **Statistical completeness**: Strand bias calculation for all variants
