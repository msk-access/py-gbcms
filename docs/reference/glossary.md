# Glossary

Technical terms used throughout the documentation.

## Metrics

| Term | Definition |
|:-----|:-----------|
| **DP** | Total depth — all mapped, quality-filtered reads overlapping the variant anchor position, including REF, ALT, and 'neither'. `DP ≥ RD + AD`. |
| **VAF** | Variant Allele Frequency — `AD / (RD + AD)` |
| **Strand Bias** | Fisher's exact test for read direction imbalance |
| **AD** | Alternate allele depth (supporting reads) |
| **RD** | Reference allele depth |

## File Formats

| Term | Definition |
|:-----|:-----------|
| **VCF** | Variant Call Format — standard variant file |
| **MAF** | Mutation Annotation Format — annotation-rich variant file |
| **BAM** | Binary Alignment Map — compressed alignment file |
| **BAI** | BAM Index — enables random access to BAM |
| **FASTA** | Reference genome sequence file |
| **FAI** | FASTA Index — enables random access to FASTA |

## Quality Scores

| Term | Definition |
|:-----|:-----------|
| **MAPQ** | Mapping Quality — confidence in read alignment |
| **BASEQ** | Base Quality — confidence in base call |

## cfDNA Terms

| Term | Definition |
|:-----|:-----------|
| **cfDNA** | Cell-free DNA — circulating DNA in plasma |
| **ctDNA** | Circulating tumor DNA — tumor-derived cfDNA |
| **Duplex** | Reads from both strands of original molecule |

## Fragment Metrics

| Term | Definition |
|:-----|:-----------|
| **DPF** | Fragment depth — total unique fragments (including discarded) |
| **RDF** | Reference fragment count — fragments resolved to REF |
| **ADF** | Alternate fragment count — fragments resolved to ALT |
| **Fragment consensus** | Quality-weighted method to resolve R1/R2 disagreements within a fragment |

## Validation Status

| Status | Meaning |
|:-------|:--------|
| **PASS** | Variant validated against reference, ready for counting |
| **PASS_WARN_REF_CORRECTED** | REF allele ≥90% match; corrected to FASTA REF |
| **PASS_WARN_HOMOPOLYMER_DECOMP** | Variant spans a homopolymer; dual-counted with corrected allele (corrected won) |
| **PASS_MULTI_ALLELIC** | Variant overlaps another variant at the same locus; sibling ALT exclusion active |
| **REF_MISMATCH** | REF allele does not match the reference genome at the stated position |
| **FETCH_FAILED** | Could not fetch the reference region (chromosome not found, etc.) |

## Related

- [Architecture](architecture.md) — System design
- [Input Formats](input-formats.md) — File specifications
