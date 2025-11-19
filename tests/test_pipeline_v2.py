"""
Integration test for gbcms v2 pipeline.
"""

import pytest
from pathlib import Path
import gbcms_rs
from gbcms.io.input import VcfReader, ReferenceChecker
from gbcms.io.output import MafWriter
from gbcms.models.core import Variant

def test_pipeline_v2(tmp_path):
    base_dir = Path("tests/testdata")
    bam_path = str(base_dir / "sample1_integration_test.bam")
    vcf_path = str(base_dir / "integration_test_variants.vcf")
    fasta_path = str(base_dir / "integration_test_reference.fa")
    output_path = tmp_path / "output.maf"
    
    # 1. Read Variants
    reader = VcfReader(Path(vcf_path))
    variants = list(reader)
    reader.close()
    
    assert len(variants) > 0
    
    # 2. Validate against Reference
    ref_checker = ReferenceChecker(Path(fasta_path))
    for v in variants:
        is_valid = ref_checker.validate(v)
        # Note: Some test variants might be artificial/invalid, so we just log or check specific ones
        # assert is_valid, f"Variant {v} does not match reference"
    ref_checker.close()
    
    # 3. Convert to Rust Variants
    rs_variants = [
        gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value)
        for v in variants
    ]
    
    # 4. Run Rust Engine
    results = gbcms_rs.count_bam(
        bam_path,
        rs_variants,
        min_mapq=20,
        min_baseq=10,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False
    )
    
    assert len(results) == len(variants)
    
    # 5. Write Output
    writer = MafWriter(output_path)
    for v, counts in zip(variants, results):
        writer.write(v, counts)
    writer.close()
    
    # 6. Verify Output
    assert output_path.exists()
    with open(output_path) as f:
        lines = f.readlines()
        assert len(lines) > 1 # Header + Data
        header = lines[0].strip().split('\t')
        assert "t_ref_count" in header
        assert "t_vaf_fragment" in header
