"""
Tests for MafReader after normalization refactor.

MafReader now yields raw MAF coordinates — it no longer performs anchor base
resolution (that's handled by Rust prepare_variants()).  These tests verify
that the reader correctly parses MAF rows and converts them to internal
coordinates via CoordinateKernel.maf_to_internal().
"""

import tempfile
import unittest
from pathlib import Path

from gbcms.io.input import MafReader
from gbcms.models.core import VariantType


class TestMafReader(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_path = Path(self.test_dir.name)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_maf_snp(self):
        """SNP: raw MAF coords pass through directly."""
        maf_path = self.base_path / "test_snp.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\n"
            )
            f.write("Gene1\tchr1\t100\t100\tA\tT\n")

        reader = MafReader(maf_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]

        self.assertEqual(v.chrom, "1")  # chr prefix stripped
        self.assertEqual(v.pos, 99)  # 1-based → 0-based
        self.assertEqual(v.ref, "A")
        self.assertEqual(v.alt, "T")
        self.assertEqual(v.variant_type, VariantType.SNP)

    def test_maf_insertion_raw_coords(self):
        """Insertion: MafReader yields raw MAF coords (Ref='-').

        Anchor base resolution is now done by Rust prepare_variants().
        """
        maf_path = self.base_path / "test_ins.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\n"
            )
            f.write("Gene1\tchr1\t5\t6\t-\tG\n")

        reader = MafReader(maf_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]

        self.assertEqual(v.chrom, "1")
        # maf_to_internal converts MAF insertion: pos = start_pos - 1 (0-based)
        self.assertEqual(v.ref, "-")
        self.assertEqual(v.alt, "G")
        self.assertEqual(v.variant_type, VariantType.INSERTION)

    def test_maf_deletion_raw_coords(self):
        """Deletion: MafReader yields raw MAF coords (Alt='-').

        Anchor base resolution is now done by Rust prepare_variants().
        """
        maf_path = self.base_path / "test_del.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\n"
            )
            f.write("Gene1\tchr1\t6\t6\tT\t-\n")

        reader = MafReader(maf_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]

        self.assertEqual(v.chrom, "1")
        self.assertEqual(v.ref, "T")
        self.assertEqual(v.alt, "-")
        self.assertEqual(v.variant_type, VariantType.DELETION)

    def test_maf_complex_raw_coords(self):
        """Complex variant: raw coords, no anchor resolution."""
        maf_path = self.base_path / "test_complex.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\n"
            )
            f.write("Gene1\tchr1\t6\t7\tCG\tTAA\n")

        reader = MafReader(maf_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]

        self.assertEqual(v.chrom, "1")
        self.assertEqual(v.ref, "CG")
        self.assertEqual(v.alt, "TAA")
        self.assertEqual(v.variant_type, VariantType.COMPLEX)

    def test_maf_metadata_preserved(self):
        """Verify that original MAF columns are preserved in metadata."""
        maf_path = self.base_path / "test_meta.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
            )
            f.write("TP53\tchr17\t7577120\t7577120\tC\tT\tSAMPLE_001\n")

        reader = MafReader(maf_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]

        # Metadata should contain all original columns
        self.assertIn("Hugo_Symbol", v.metadata)
        self.assertEqual(v.metadata["Hugo_Symbol"], "TP53")
        self.assertEqual(v.metadata["Tumor_Sample_Barcode"], "SAMPLE_001")


if __name__ == "__main__":
    unittest.main()
