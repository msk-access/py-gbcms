"""
Benchmark script for gbcms v2 (Rust) vs v1 (Python).
"""

import sys
import time
from pathlib import Path

# Add src to path to import gbcms_rs
sys.path.append("src")

from gbcms.config import Config
from gbcms.processor import VariantProcessor
from gbcms_v2.models.core import Variant, VariantType

import gbcms_rs


def benchmark_rust(bam_path: str, variants: list[Variant]):
    print(f"Benchmarking Rust engine with {len(variants)} variants...")
    start = time.time()

    # Convert variants to Rust-compatible objects
    rs_variants = [
        gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants
    ]

    results = gbcms_rs.count_bam(
        bam_path,
        rs_variants,
        min_mapq=20,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
    )

    end = time.time()
    print(f"Rust engine finished in {end - start:.4f} seconds")
    return results


def benchmark_python(bam_path: str, vcf_path: str, fasta_path: str):
    print("Benchmarking Python engine...")
    start = time.time()

    config = Config(
        fasta_file=fasta_path,
        bam_files={"sample": bam_path},
        variant_files=[vcf_path],
        output_file="/dev/null",
        input_is_vcf=True,
    )

    processor = VariantProcessor(config)
    processor.process()

    end = time.time()
    print(f"Python engine finished in {end - start:.4f} seconds")


if __name__ == "__main__":
    base_dir = Path("tests/testdata")
    bam_path = str(base_dir / "sample1_integration_test.bam")
    vcf_path = str(base_dir / "integration_test_variants.vcf")
    fasta_path = str(base_dir / "integration_test_reference.fa")

    print(f"Using test data from: {base_dir}")

    # Load variants using legacy loader for fairness/convenience
    from gbcms.variant import VariantLoader

    loader = VariantLoader()
    legacy_variants = loader.load_vcf(vcf_path)

    # Convert to new model
    new_variants = []

    for v in legacy_variants:
        # Legacy variants are already 0-based for start?
        # VariantEntry: pos is 0-based.
        # We need to map types.
        v_type = VariantType.SNP
        if v.insertion:
            v_type = VariantType.INSERTION
        elif v.deletion:
            v_type = VariantType.DELETION

        new_variants.append(
            Variant(chrom=v.chrom, pos=v.pos, ref=v.ref, alt=v.alt, variant_type=v_type)
        )

    # Run benchmarks
    benchmark_rust(bam_path, new_variants)
    benchmark_python(bam_path, vcf_path, fasta_path)
