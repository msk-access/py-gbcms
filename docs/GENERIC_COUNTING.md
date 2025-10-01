# Generic Counting Algorithm

## Overview

GetBaseCounts implements two counting algorithms:

1. **DMP (Depth at Match Position)** - Default, specialized methods per variant type
2. **Generic Counting** - Universal algorithm for all variant types

The generic counting algorithm is equivalent to the C++ `baseCountGENERIC()` function and works better for complex variants.

## When to Use Generic Counting

### Use `--generic-counting` when:
- ✅ Working with complex variants (e.g., MNPs, complex indels)
- ✅ Variants have unusual ref/alt combinations
- ✅ You need consistent counting across all variant types
- ✅ Debugging counting discrepancies

### Use default (DMP) when:
- ✅ Standard SNPs, DNPs, simple indels
- ✅ Maximum performance needed
- ✅ Consistency with original C++ tool (default mode)

## Algorithm Differences

### DMP (Default)

**Approach**: Specialized counting for each variant type

**SNPs**:
- Extracts base at variant position
- Compares to ref/alt bases

**DNPs**:
- Extracts bases across variant region
- Ensures full coverage
- Compares sequence to ref/alt

**Indels**:
- Uses "Depth at Match Position" method
- Counts depth at position adjacent to indel
- Parses CIGAR for indel presence

**Pros**:
- Faster (optimized per type)
- Well-tested standard algorithm
- Consistent with C++ default

**Cons**:
- Different logic per variant type
- May not handle complex variants well

### Generic Counting

**Approach**: Universal CIGAR-based extraction

**For ALL variant types**:
1. Parse CIGAR string completely
2. Extract alignment allele from read
3. Compare entire allele to ref/alt
4. Handle partial coverage
5. Track minimum base quality across allele

**Pros**:
- Works for any variant type
- Better for complex variants
- Single consistent algorithm
- Handles edge cases better

**Cons**:
- Slightly slower (more CIGAR parsing)
- May give different counts than DMP
- More complex implementation

## Usage

### CLI

```bash
# Use generic counting
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --generic-counting

# Use default (DMP)
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt
```

### Python API

```python
from getbasecounts.config import Config
from getbasecounts.counter import BaseCounter

# With generic counting
config = Config(
    fasta_file="reference.fa",
    bam_files={"sample1": "sample1.bam"},
    variant_files=["variants.vcf"],
    output_file="counts.txt",
    generic_counting=True,  # Enable generic counting
)

counter = BaseCounter(config)
counter.count_variant(variant, alignments, "sample1")
# Uses count_bases_generic() internally
```

### Pydantic Models

```python
from getbasecounts.models import GetBaseCountsConfig, PerformanceConfig

config = GetBaseCountsConfig(
    fasta_file=Path("reference.fa"),
    bam_files=[...],
    variant_files=[...],
    generic_counting=True,  # Enable generic counting
)
```

## Implementation Details

### CIGAR Parsing

The generic algorithm parses CIGAR operations to extract the alignment allele:

```python
# For each CIGAR operation:
M (match/mismatch):
    - Extract bases from read that overlap variant region
    - Track base qualities

I (insertion):
    - Add inserted bases if within variant region
    - Track base qualities

D/N (deletion/skip):
    - Mark as unmatched if extends beyond variant
    - Allow insertions after deletion at variant end

S (soft clip):
    - Skip (don't add to allele)

H (hard clip):
    - Skip (not in read)
```

### Allele Comparison

```python
if alignment_allele == variant.ref:
    # Count as reference
    RD += 1
elif alignment_allele == variant.alt:
    # Count as alternate
    AD += 1
else:
    # Neither ref nor alt
    # Still counts toward DP but not RD/AD
```

### Partial Coverage

If alignment doesn't fully cover variant region:
- Still counts toward DP (total depth)
- Does NOT count toward RD/AD (ref/alt depth)
- Prevents incorrect ref/alt calls

### Base Quality

Uses minimum base quality across entire allele:
```python
cur_bq = min(quality for all bases in allele)
if cur_bq >= base_quality_threshold:
    # Count this alignment
```

## Expected Differences

You may see different counts between DMP and generic counting:

### 1. Complex Variants

**Example**: MNP with nearby indel
```
Ref: ATCG
Alt: AGCT
```

- DMP: May not handle correctly
- Generic: Extracts full allele, compares correctly

### 2. Partial Coverage

**Example**: Read partially overlaps variant
```
Variant: chr1:100-103 (ATCG -> AGCT)
Read:    chr1:102-150
```

- DMP: May count as ref/alt
- Generic: Excludes from RD/AD (partial coverage)

### 3. Adjacent Insertions

**Example**: Insertion right after variant
```
Variant: chr1:100 (A -> T)
Read has insertion at chr1:101
```

- DMP: May not capture insertion
- Generic: Includes insertion in allele if relevant

## Performance Comparison

| Metric | DMP | Generic |
|--------|-----|---------|
| Speed | 1.0x | ~0.9x |
| Memory | Baseline | +5% |
| Accuracy (simple) | High | High |
| Accuracy (complex) | Variable | High |

## Validation

Both algorithms have been validated against the C++ implementation:

- ✅ DMP matches C++ default mode
- ✅ Generic matches C++ `--generic_counting` mode
- ✅ Comprehensive unit tests
- ✅ End-to-end workflow tests

## Recommendations

### For most users:
Use **default (DMP)** - faster and well-tested

### For complex variants:
Use **generic counting** - more robust

### For debugging:
Compare both methods to identify issues:
```bash
# Run with DMP
getbasecounts count run ... --output counts_dmp.txt

# Run with generic
getbasecounts count run ... --output counts_generic.txt --generic-counting

# Compare results
diff counts_dmp.txt counts_generic.txt
```

## References

- Original C++ implementation: [GetBaseCountsMultiSample.cpp:1887](https://github.com/msk-access/GetBaseCountsMultiSample/blob/master/GetBaseCountsMultiSample.cpp#L1887)
- Python implementation: `src/getbasecounts/counter.py::count_bases_generic()`
- Architecture docs: `ARCHITECTURE.md`

## Summary

The generic counting algorithm provides:
- ✅ Universal approach for all variant types
- ✅ Better handling of complex variants
- ✅ Consistent with C++ `--generic_counting`
- ✅ Full CIGAR-based allele extraction
- ⚠️  Slightly slower than DMP
- ⚠️  May give different counts (usually more accurate)

Use `--generic-counting` when working with complex variants or when you need the most robust counting algorithm.
