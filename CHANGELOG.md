# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-11-20

### Added
- **Rust Backend**: Complete rewrite of the counting engine in Rust for significant performance improvements.
- **Strand-Specific Counts**: Output now includes forward/reverse counts for reads and fragments by default.
- **VAF/FAF**: Variant Allele Fraction (VAF) and Fragment Allele Fraction (FAF) are now calculated and included in the output.
- **New Filters**: Added `filter_qc_failed`, `filter_improper_pair`, and `filter_indel` options.
- **CLI Improvements**:
    - Direct `gbcms run` command.
    - Explicit sample ID support via `--bam sample_id:path`.
    - Output suffix customization via `--suffix`.
- **MAF Column Preservation**: Input MAF columns are now preserved in the output.

### Changed
- **Architecture**: Moved from pure Python to a hybrid Python/Rust architecture.
- **Output Format**:
    - VCF `FORMAT` fields updated to include strand counts (e.g., `RD` is now `ref_fwd,ref_rev`).
    - MAF output now includes standard `t_ref_count_forward`, `t_ref_count_reverse`, etc.
- **CLI**: Deprecated `python -m gbcms.cli` in favor of `gbcms`.

### Removed
- **Legacy Python Counting**: The old pure Python counting implementation has been removed.
