# gbcms run

Count alleles at variant positions across one or more BAM files.

## Synopsis

```bash
gbcms run [OPTIONS] --variants <FILE> --bam <NAME:PATH>... --fasta <FILE>
```

## Required Arguments

| Option | Description |
|:-------|:------------|
| `--variants`, `-v` | [VCF or MAF](../reference/input-formats.md) file with variant positions (`.vcf`, `.vcf.gz`, `.vcf.bgz`, or `.maf`). Unsupported extensions are rejected immediately. |
| `--bam`, `-b` | BAM file path (can repeat). Optionally prefix with `name:` for sample naming, e.g. `--bam tumor:tumor.bam`. If no name given, the filename stem is used. |
| `--bam-list`, `-L` | File containing BAM paths (one per line, optionally `sample_name path`). Alternative to repeated `--bam`. |
| `--fasta`, `-f` | Reference FASTA file (with .fai index) |
| `--lenient-bam` | Skip missing `--bam` paths and continue with remaining samples (default: exit immediately on first missing BAM). Note: a missing `--bam-list` file always fails regardless. |

## Output Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--output-dir`, `-o` | *required* | Output directory |
| `--format` | `vcf` | Output format (`vcf` or `maf`) |
| `--suffix` | `''` | Suffix for output filenames |
| `--column-prefix` | `''` | Prefix for gbcms count columns in MAF output. Only letters, digits, and underscores (`[A-Za-z0-9_]`) are allowed; invalid characters exit immediately. |
| `--preserve-barcode` | `false` | Keep original Tumor_Sample_Barcode from input MAF. No-op (with warning) when input is not MAF. |
| `--show-normalization` | `false` | Append `norm_*` columns showing left-aligned coordinates |
| `--context-padding` | `5` | Minimum flanking bases for haplotype construction. Range **1–50**, enforced at parse time. Auto-increased in repeat regions when `--adaptive-context` is enabled. |
| `--adaptive-context` | `true` | Dynamically increase context padding in [tandem repeat regions](../reference/variant-normalization.md#adaptive-context-padding) |
| `--threads` | `1` | Number of threads |

## mFSD Options

Mutant Fragment Size Distribution (mFSD) analysis compares insert-size distributions
for REF- vs ALT-classified fragments at each variant position, enabling detection of
short-fragment enrichment associated with tumor-derived cfDNA
(see [mFSD Metrics](../reference/counting-metrics.md#mfsd)).

| Option | Default | Description |
|:-------|:--------|:------------|
| `--mfsd` | `false` | Enable mFSD analysis. Adds 34 mFSD columns (KS test, LLR, mean sizes, pairwise comparisons, derived metrics) to MAF output and 7 `MFSD_*` INFO fields to VCF. |
| `--mfsd-parquet` | `false` | Write a companion `<sample>.fsd.parquet` with per-variant raw fragment size arrays (`ref_sizes`, `alt_sizes`). Enables downstream visualizations. **Requires `--mfsd`**. |

!!! tip
    To generate both summary statistics and raw Parquet data in one run:
    ```bash
    gbcms run --mfsd --mfsd-parquet --format maf \\
        --variants variants.maf --bam tumor:tumor.bam --fasta hg19.fa -o ./results
    ```
    This produces `<sample>.maf` (with 34 mFSD columns) and `<sample>.fsd.parquet`
    (raw fragment sizes for visualization).

## Filtering Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--min-mapq` | `20` | Minimum MAPQ |
| `--min-baseq` | `20` | Minimum BASEQ |
| `--filter-duplicates` | `true` | Filter duplicate reads |
| `--filter-secondary` | `false` | Filter secondary alignments |
| `--filter-supplementary` | `false` | Filter supplementary alignments |
| `--filter-qc-failed` | `false` | Filter QC failed reads |
| `--filter-improper-pair` | `false` | Filter improperly paired reads |
| `--filter-indel` | `false` | Filter reads with indels |
| `--fragment-qual-threshold` | `10` | Quality difference threshold for fragment consensus (see [Fragment Counting](../reference/counting-metrics.md#fragment-counting)) |

## Debugging Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `--verbose`, `-V` | `false` | Enable verbose debug logging |
| `--trace`, `-T` | `false` | Enable per-read Rust trace logging (slow). Implies `--verbose`. Shows detailed per-read classification diagnostics. |

## Advanced: Alignment Backend

By default, Phase 3 allele classification uses Smith-Waterman (SW). An alternative PairHMM backend is available for probabilistic alignment:

| Option | Default | Description |
|:-------|:--------|:------------|
| `--alignment-backend` | `sw` | Alignment backend for Phase 3: `sw` (Smith-Waterman) or `hmm` (PairHMM). Invalid values are rejected at parse time. |
| `--llr-threshold` | `2.3` | PairHMM log-likelihood ratio threshold for confident calls (≈ ln(10)) |
| `--gap-open-prob` | `1e-4` | PairHMM gap-open probability for non-repeat regions |
| `--gap-extend-prob` | `0.1` | PairHMM gap-extend probability for non-repeat regions |
| `--repeat-gap-open-prob` | `1e-2` | PairHMM gap-open probability for tandem repeat regions |
| `--repeat-gap-extend-prob` | `0.5` | PairHMM gap-extend probability for tandem repeat regions |

!!! tip "When to use PairHMM"
    The PairHMM backend uses base quality probabilities directly in alignment scoring, making it more sensitive in low-quality or noisy regions. For most workflows, the default SW backend is recommended. Use `--alignment-backend hmm` when you need probabilistic confidence scores (LLR) instead of edit-distance margins.

## Examples

### Single BAM

```bash
gbcms run \
    --variants mutations.vcf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --output-dir results/
```

### Multiple BAMs

```bash
gbcms run \
    --variants mutations.maf \
    --bam tumor:tumor.bam \
    --bam normal:normal.bam \
    --fasta reference.fa \
    --format maf
```

### With Filtering

```bash
gbcms run \
    --variants mutations.vcf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --filter-duplicates \
    --min-mapq 30
```

### With Normalization Columns

```bash
gbcms run \
    --variants mutations.maf \
    --bam sample:sample.bam \
    --fasta reference.fa \
    --format maf \
    --show-normalization
```

## Output

The command produces a VCF or MAF file with:

- **Allele counts** (reference and alternate depth)
- **VAF** (variant allele frequency)
- **Strand bias** (Fisher's exact test)
- **Fragment counts** (deduplicated)
- **Validation status** (`PASS`, `PASS_WARN_HOMOPOLYMER_DECOMP`, `REF_MISMATCH`, or `FETCH_FAILED`)
- **Normalization columns** (with `--show-normalization`): left-aligned position, REF, and ALT

## Related

- [Quick Start](../getting-started/quickstart.md) — Common patterns
- [gbcms normalize](normalize.md) — Standalone normalization (no counting)
- [Nextflow Pipeline](../nextflow/index.md) — For many samples
- [Input Formats](../reference/input-formats.md) — VCF/MAF specs
- [Variant Counting](../reference/allele-classification.md) — How each variant type is counted
