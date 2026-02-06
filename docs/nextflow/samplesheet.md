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

## Related

- [Parameters](parameters.md) — All pipeline options
- [Examples](examples.md) — Common usage patterns
- [Troubleshooting](../resources/troubleshooting.md) — Common issues
