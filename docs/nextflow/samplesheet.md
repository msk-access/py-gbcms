# Samplesheet Format

The samplesheet is a CSV file defining input samples for the Nextflow pipeline.

## Basic Format

```csv
sample,bam,bai
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
sample2,/path/to/sample2.bam,
sample3,/path/to/sample3.bam,/path/to/sample3.bam.bai
```

## Columns

| Column | Required | Description |
|:-------|:---------|:------------|
| `sample` | Yes | Sample identifier (used in output filenames) |
| `bam` | Yes | Path to BAM file |
| `bai` | No | Path to BAM index (auto-discovers `.bam.bai` or `.bai`) |
| `suffix` | No | Per-sample output suffix |
| `tsb` | No | Tumor_Sample_Barcode pattern(s) for [MAF filtering](#multi-sample-maf-filtering) |

## BAI Auto-Discovery

If `bai` column is empty, the pipeline checks for:

1. `<bam>.bai` (e.g., `sample.bam.bai`)
2. `<bam_without_extension>.bai` (e.g., `sample.bai`)

!!! tip
    Leave the `bai` column empty to use auto-discovery.

## Per-Sample Suffix

For samples with multiple BAM types:

```csv
sample,bam,bai,suffix
sample1,/path/to/sample1.duplex.bam,,-duplex
sample1,/path/to/sample1.simplex.bam,,-simplex
sample1,/path/to/sample1.unfiltered.bam,,-unfiltered
```

Output files: `sample1-duplex.maf`, `sample1-simplex.maf`, `sample1-unfiltered.maf`

## Multi-Sample MAF Filtering

When using `--filter_by_sample` with a multi-sample MAF as `--variants`, each sample's variants are filtered by `Tumor_Sample_Barcode`.

### Exact Match (default)

When no `tsb` column is provided, the `sample` column is matched **exactly** against `Tumor_Sample_Barcode`:

```csv
sample,bam
P-0012345-T01-IM7,/path/to/tumor.bam
P-0067890-T01-IM7,/path/to/tumor2.bam
```

### Regex Match

Use the `tsb` column for pattern matching (e.g., patient-level filtering):

```csv
sample,bam,tsb
patient_A,/path/to/A.bam,P-0012345
patient_B,/path/to/B.bam,P-0067890
```

`P-0012345` matches any `Tumor_Sample_Barcode` containing that substring (e.g., `P-0012345-T01-IM7`, `P-0012345-T02-IM7`).

### Comma-Separated Multi-Select

Select specific samples with comma-separated patterns:

```csv
sample,bam,tsb
tumor_pair,/path/to/A.bam,"P-0012345-T01-IM7,P-0012345-T02-IM7"
```

Duplicate rows (matched by multiple patterns) are automatically deduplicated.

!!! note
    Samples with 0 matching variants are skipped and reported in `pipeline_summary.tsv`.

## Related

- [Parameters](parameters.md) — All pipeline options
- [Examples](examples.md) — Common usage patterns
- [Troubleshooting](../resources/troubleshooting.md) — Common issues
