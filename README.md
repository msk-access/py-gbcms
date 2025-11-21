# gbcms: GetBaseCounts Multi-Sample for variant counting from BAM files

An adaptation of the original gbcms tool to be used as a multi-sample variant counting system.

Original gbcms: https://github.com/zengzheng123/GetBaseCountsMultiSample

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

## 📚 Documentation

**Full documentation is available at: https://msk-access.github.io/py-gbcms/**

Quick Links:
- [Quick Start Guide](https://msk-access.github.io/py-gbcms/quick-start.html)
- [Installation Guide](https://msk-access.github.io/py-gbcms/INSTALLATION.html)
- [CLI Reference](https://msk-access.github.io/py-gbcms/CLI_FEATURES.html)
- [Advanced Usage](https://msk-access.github.io/py-gbcms/advanced-usage.html)
- [Architecture](https://msk-access.github.io/py-gbcms/architecture.html)

**Contributing to Documentation:**
Documentation lives in the [`gh-pages` branch](https://github.com/msk-access/py-gbcms/tree/gh-pages). 
See the [gh-pages README](https://github.com/msk-access/py-gbcms/blob/gh-pages/README.md) for contribution guidelines.

## Installation

For detailed installation instructions, see the [Installation Guide](https://msk-access.github.io/py-gbcms/INSTALLATION.html).

### Quick Install (from source)

```bash
# Requires Rust toolchain
pip install .
```

### Docker

```bash
docker pull msk-access/gbcms:latest
# or build locally
docker build -t gbcms:latest .
```

## Quick Start


### Count variants from VCF with BAM alignments
```bash
gbcms run --fasta reference.fa --bam sample1.bam --variants variants.vcf --output-dir results/
```

### Count variants from MAF with BAM alignments (sample-agnostic MAF output)
```bash
gbcms run --fasta reference.fa --bam sample1.bam --variants variants.maf --output-dir results/ --format maf
```
### Process multiple BAMs using a File of Files (FoF)
```bash
gbcms run --fasta reference.fa --bam-list bam_list.txt --variants variants.vcf --output-dir results/
```
## Advanced Usage

### Explicit Sample IDs
You can explicitly specify sample IDs for your BAM files. This ID will be used in the output filename and internal file headers.

**CLI Argument:**

```bash
gbcms run ... --bam "MySampleID:/path/to/sample.bam"
```

**BAM List File:**
The BAM list file supports a two-column format (whitespace separated):
```text
SampleID1   /path/to/sample1.bam
SampleID2   /path/to/sample2.bam
```

### Output Suffix
You can append a custom suffix to the output filenames using the `--suffix` flag.
```bash
gbcms run ... --suffix .genotyped
# Output: {SampleID}.genotyped.vcf
```

## Output Format Details

gbcms provides detailed strand-specific counts and allele fractions by default.

### VCF Output
The following `FORMAT` fields are added:

- **DP**: Total Depth (`ref_total,alt_total`)
- **RD**: Reference Read Depth (`fwd,rev`)
- **AD**: Alternate Read Depth (`fwd,rev`)
- **RDF**: Reference Fragment Depth (`fwd,rev`)
- **ADF**: Alternate Fragment Depth (`fwd,rev`)
- **VAF**: Variant Allele Fraction (Read Level)
- **FAF**: Variant Allele Fraction (Fragment Level)

### MAF Output
The following columns are added to the MAF output:

- `t_ref_count_forward`, `t_ref_count_reverse`
- `t_alt_count_forward`, `t_alt_count_reverse`
- `t_ref_count_fragment_forward`, `t_ref_count_fragment_reverse`
- `t_alt_count_fragment_forward`, `t_alt_count_fragment_reverse`
- `t_vaf` (Read Level VAF)
- `t_vaf_fragment` (Fragment Level VAF)

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
