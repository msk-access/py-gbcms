# Input Formats

py-gbcms accepts VCF and MAF files as variant input.

## VCF (Variant Call Format)

Standard VCF format with required fields:

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    12345   .       A       T       .       PASS    .
chr2    67890   .       G       C       .       PASS    .
```

### Requirements

- Tab-separated
- `#CHROM`, `POS`, `REF`, `ALT` columns required
- 1-based positions

## MAF (Mutation Annotation Format)

Standard MAF format with required columns:

```
Hugo_Symbol  Chromosome  Start_Position  End_Position  Reference_Allele  Tumor_Seq_Allele2
TP53         chr17       7577120         7577120       C                 T
KRAS         chr12       25398284        25398284      G                 A
```

### Required Columns

| Column | Description |
|:-------|:------------|
| `Chromosome` | Chromosome name |
| `Start_Position` | 1-based start position |
| `Reference_Allele` | Reference allele |
| `Tumor_Seq_Allele2` | Alternate allele |

## Reference FASTA

- Must have corresponding `.fai` index
- Chromosome names must match VCF/MAF

```bash
# Create index if missing
samtools faidx reference.fa
```

## BAM Requirements

- Must have corresponding `.bai` index
- Coordinate-sorted
- Chromosome names must match reference

## Related

- [CLI Run Command](../cli/run.md) — Usage examples
- [Architecture](architecture.md) — How counting works
- [Glossary](glossary.md) — Term definitions
