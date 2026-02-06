# CLI Reference

The `gbcms` command-line interface provides a single powerful command for variant counting.

## Commands

| Command | Description |
|:--------|:------------|
| [**run**](run.md) | Count alleles at variant positions |

## Quick Example

```bash
gbcms run \
    --variants mutations.maf \
    --bam sample1:sample1.bam \
    --bam sample2:sample2.bam \
    --fasta reference.fa \
    --output-dir results/
```

## Getting Help

```bash
gbcms --help
gbcms run --help
```

## Related

- [Quick Start](../getting-started/quickstart.md) — Common usage patterns
- [Nextflow Pipeline](../nextflow/index.md) — For processing many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specifications
