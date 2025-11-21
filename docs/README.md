# Introduction

**gbcms** (Get Base Counts Multi-Sample) is a high-performance tool for extracting base counts and variant metrics from BAM files. It is designed for accuracy, speed, and ease of use in bioinformatics pipelines.

## Key Features

*   **Rust-Powered Engine**: Core counting logic is implemented in Rust for maximum performance and memory efficiency.
*   **Accurate Variant Handling**: 
    *   Full support for VCF and MAF input formats.
    *   Rigorous normalization of indels and complex variants using the reference genome.
    *   Validation of input variants against the reference FASTA.
*   **Comprehensive Metrics**:
    *   Standard depth (DP), reference depth (RD), and alternate depth (AD).
    *   Strand-specific counts (Forward/Reverse).
    *   Fragment-level counts (deduplicated by fragment ID).
    *   Variant Allele Fractions (VAF) at read and fragment levels.
    *   **Strand Bias Statistics**: Fisher's Exact Test p-values (`FSB_PVAL`) and Odds Ratios (`FSB_OR`) calculated at the fragment level.
*   **Modern CLI**: User-friendly command-line interface with rich output and progress bars.
*   **Flexible Input**: Support for single BAMs, lists of BAMs (File of Files), and automatic sample naming.

## Architecture

gbcms uses a hybrid Python/Rust architecture:
*   **Python (v2)**: Handles CLI, input parsing (VCF/MAF), orchestration, and output formatting.
*   **Rust**: Performs the computationally intensive task of iterating over reads, piling up bases, and calculating statistics.

## Next Steps

*   [Install gbcms](INSTALLATION.md)
*   [Learn how to use the CLI](CLI_FEATURES.md)
*   [Understand Input/Output formats](INPUT_OUTPUT.md)
