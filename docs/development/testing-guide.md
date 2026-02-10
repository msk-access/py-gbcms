# Testing Guide

This guide covers running tests, adding new tests, and accuracy validation for py-gbcms.

## Running Tests

### Quick Test

```bash
# Run all tests
pytest -v

# Run with coverage
pytest --cov=gbcms --cov-report=html

# Run specific test file
pytest tests/test_accuracy.py -v
```

### Test Categories

| Category | Files | Purpose |
|:---------|:------|:--------|
| Accuracy | `test_accuracy.py` | SNP, indel, complex variant counting |
| Shifted Indels | `test_shifted_indels.py` | Windowed indel detection (±5bp), 3-layer safeguards |
| Complex Masking | `test_fuzzy_complex.py` | Quality-aware masked comparison, ambiguity detection |
| CLI | `test_cli_sample_id.py` | Command-line parsing |
| Filters | `test_filters.py` | Read filtering logic |
| MAF | `test_maf_*.py` | MAF column preservation |
| Pipeline | `test_pipeline_v2.py` | End-to-end workflow |
| Strand | `test_strand_counts.py` | Strand-specific counts |

---

## Test Structure

```
tests/
├── test_accuracy.py           # Variant type accuracy
├── test_cli_sample_id.py      # CLI argument parsing
├── test_filters.py            # Read filtering
├── test_fuzzy_complex.py      # Quality-aware masked complex matching
├── test_maf_preservation.py
├── test_maf_reader.py
├── test_pipeline_v2.py
├── test_shifted_indels.py     # Windowed indel detection (±5bp)
└── test_strand_counts.py
```

---

## Writing Tests

### Basic Test Template

```python
import pytest
from pathlib import Path

def test_my_feature(tmp_path):
    """Test description."""
    # Arrange
    input_file = tmp_path / "input.txt"
    input_file.write_text("test data")
    
    # Act
    result = my_function(input_file)
    
    # Assert
    assert result.success
    assert result.count == 42
```

### Accuracy Test Template

```python
def test_snp_accuracy():
    """Verify SNP counting against known BAM."""
    # Create variant
    variant = Variant("chr1", 100, "A", "T", "SNP")
    
    # Run counting
    results = count_bam(bam_path, [variant], fasta_path)
    
    # Validate
    assert results[0].ref_count == 50
    assert results[0].alt_count == 10
```

---

## Manual Validation

### Using samtools for Spot-Check

```bash
# Check counts at specific position
samtools mpileup -r chr1:100-100 -q 20 \
    -f ref.fa sample.bam 2>/dev/null | \
    awk '{print "DP="$4}'
```

### Comparing with gbcms Output

```bash
# Run gbcms
gbcms run -v variants.maf -b sample.bam -f ref.fa -o output/

# Check output
awk -F'\t' 'NR==2 {print "REF="$41, "ALT="$42}' output/*.maf
```

---

## Accuracy Validation

### Variant Types Tested

| Type | Test | Status |
|:-----|:-----|:------:|
| SNP | `test_snp_accuracy` | ✅ |
| Insertion | `test_insertion_accuracy` | ✅ |
| Deletion | `test_deletion_accuracy` | ✅ |
| Complex | `test_complex_accuracy` | ✅ |
| MNP | `test_mnp_accuracy` | ✅ |
| Shifted Indels | `test_shifted_indels.py` (15 cases) | ✅ |
| Complex Masking | `test_fuzzy_complex.py` (14 cases) | ✅ |

### Real-World Validation

```bash
# Compare gbcms vs samtools for a SNP
# Position: chr1:11168293 G>A

# gbcms output
awk -F'\t' '$5=="1" && $6=="11168293" {print "REF="$41, "ALT="$42}' output.maf

# samtools output
samtools mpileup -r 1:11168293-11168293 -q 20 -f ref.fa sample.bam | \
    awk '{gsub(/\^.|\$/,"",$5); print "DP="$4, "Pileup="$5}'
```

---

## Coverage Targets

| Module | Target | Current |
|:-------|:------:|:-------:|
| cli.py | 80% | 79% |
| pipeline.py | 70% | 29% |
| io/input.py | 85% | 82% |
| io/output.py | 90% | 96% |
| models/core.py | 90% | 90% |

Run coverage report:

```bash
pytest --cov=gbcms --cov-report=html
open htmlcov/index.html
```
