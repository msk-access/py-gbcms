import unittest
import tempfile
from pathlib import Path
import pysam
from gbcms.io.input import MafReader
from gbcms.models.core import VariantType

class TestMafReader(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_path = Path(self.test_dir.name)
        
        # Create dummy FASTA
        self.fasta_path = self.base_path / "ref.fa"
        with open(self.fasta_path, "w") as f:
            f.write(">chr1\n")
            f.write("ATCGATCGATCG\n") # 012345678901
            # Pos 5 (0-based) is T. Pos 6 is C.
            
        # Index FASTA (requires samtools usually, but pysam might need it)
        # pysam.FastaFile requires an index. 
        # If we can't run samtools faidx, we might fail here.
        # Let's try to use pysam to index if possible or mock it.
        # Actually, pysam.faidx can be called from python!
        pysam.faidx(str(self.fasta_path))

    def tearDown(self):
        self.test_dir.cleanup()

    def test_maf_insertion(self):
        # Insertion of 'G' after pos 5 (T).
        # MAF: Start=5 (anchor), Ref='-', Alt='G'
        # VCF: Pos=5, Ref='T', Alt='TG'
        
        maf_path = self.base_path / "test.maf"
        with open(maf_path, "w") as f:
            f.write("Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\n")
            f.write("Gene1\tchr1\t5\t6\t-\tG\n")
            
        reader = MafReader(maf_path, self.fasta_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]
        
        self.assertEqual(v.chrom, "1") # normalized
        self.assertEqual(v.pos, 4) # 0-based index of anchor (5 in 1-based is 4 in 0-based)
        # Wait. 
        # My logic: anchor_pos_1based = 5. anchor_pos_0based = 4.
        # FASTA: 0:A, 1:T, 2:C, 3:G, 4:A, 5:T.
        # So base at 4 is 'A'.
        # VCF POS=5.
        # So VCF REF='A'. ALT='AG'.
        
        self.assertEqual(v.ref, "A")
        self.assertEqual(v.alt, "AG")
        self.assertEqual(v.variant_type, VariantType.INSERTION)

    def test_maf_deletion(self):
        # Deletion of 'T' at pos 6.
        # MAF: Start=6 (first deleted), Ref='T', Alt='-'
        # Anchor is at 5.
        # FASTA: 0:A, 1:T, 2:C, 3:G, 4:A, 5:T.
        # Base at 5 is 'T'.
        # VCF POS=5.
        # VCF REF='TT' (Anchor 'T' + Deleted 'T').
        # VCF ALT='T' (Anchor).
        
        maf_path = self.base_path / "test_del.maf"
        with open(maf_path, "w") as f:
            f.write("Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\n")
            f.write("Gene1\tchr1\t6\t6\tT\t-\n")
            
        reader = MafReader(maf_path, self.fasta_path)
        variants = list(reader)
        self.assertEqual(len(variants), 1)
        v = variants[0]
        
        self.assertEqual(v.chrom, "1")
        self.assertEqual(v.pos, 4) # 0-based index of anchor (5 in 1-based is 4 in 0-based)
        # Base at 4 is 'A'.
        # VCF REF = 'AT' (Anchor 'A' + Deleted 'T').
        # VCF ALT = 'A'.
        
        self.assertEqual(v.ref, "AT")
        self.assertEqual(v.alt, "A")
        self.assertEqual(v.variant_type, VariantType.DELETION)

if __name__ == '__main__':
    unittest.main()
