import os

import pysam
import pytest

from gbcms._rs import count_bam


@pytest.fixture
def synthetic_bam(tmp_path):
    """
    Generates a synthetic BAM file with known reads for testing accuracy.
    Reference: chr1, 1000 bases of 'A'.
    """
    bam_path = tmp_path / "synthetic.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # 1. SNP Case: chr1:100 (0-based), REF=A, ALT=T
        # 5 Forward Reads supporting REF (A)
        for i in range(5):
            a = pysam.AlignedSegment()
            a.query_name = f"read_ref_fwd_{i}"
            a.query_sequence = "A" * 10
            a.flag = 0  # Forward
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigar = ((0, 10),)  # 10M
            a.query_qualities = [30] * 10
            outf.write(a)

        # 3 Reverse Reads supporting REF (A)
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name = f"read_ref_rev_{i}"
            a.query_sequence = "A" * 10
            a.flag = 16  # Reverse
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigar = ((0, 10),)
            a.query_qualities = [30] * 10
            outf.write(a)

        # 4 Forward Reads supporting ALT (T) at pos 100
        # Sequence: AAAAA T AAAA (T at index 5)
        # Ref start 95. Index 5 corresponds to 95+5 = 100.
        for i in range(4):
            a = pysam.AlignedSegment()
            a.query_name = f"read_alt_fwd_{i}"
            a.query_sequence = "AAAAATAAAA"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigar = ((0, 10),)
            a.query_qualities = [30] * 10
            outf.write(a)

        # 2 Reverse Reads supporting ALT (T)
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_alt_rev_{i}"
            a.query_sequence = "AAAAATAAAA"
            a.flag = 16
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigar = ((0, 10),)
            a.query_qualities = [30] * 10
            outf.write(a)

        # 2. Insertion Case: chr1:200, REF=A, ALT=AT (Insertion of T)
        # Anchor at 200.
        # 2 Forward Reads supporting Insertion
        # Cigar: 5M 1I 4M. Sequence length 10.
        # Pos 200 is index 4 (0-4 are 5 bases).
        # Ref start 196. 196, 197, 198, 199, 200 (Anchor).
        # Then Insertion.
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_ins_fwd_{i}"
            a.query_sequence = "AAAAATAAAA"  # 5th base is A (anchor), then T (ins), then AAAA
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 196
            a.mapping_quality = 60
            a.cigar = ((0, 5), (1, 1), (0, 4))  # 5M 1I 4M
            a.query_qualities = [30] * 10
            outf.write(a)

        # 3. Deletion Case: chr1:300, REF=AT, ALT=A (Deletion of T)
        # Anchor at 300.
        # 3 Forward Reads supporting Deletion
        # Cigar: 5M 1D 5M.
        # Ref start 296. 296..300 (5 bases). 301 is deleted.
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name = f"read_del_fwd_{i}"
            a.query_sequence = "AAAAAAAAAA"  # 10 bases
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 296
            a.mapping_quality = 60
            a.cigar = ((0, 5), (2, 1), (0, 5))  # 5M 1D 5M
            a.query_qualities = [30] * 10
            outf.write(a)

        # 4. Filter Checks
        # Low MAPQ (should be ignored)
        a = pysam.AlignedSegment()
        a.query_name = "read_low_mapq"
        a.query_sequence = "AAAAATAAAA"  # Supports SNP T
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 95
        a.mapping_quality = 5  # Low
        a.cigar = ((0, 10),)
        a.query_qualities = [30] * 10
        outf.write(a)

        # Low BaseQ (should be ignored)
        a = pysam.AlignedSegment()
        a.query_name = "read_low_baseq"
        a.query_sequence = "AAAAATAAAA"  # Supports SNP T
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 95
        a.mapping_quality = 60
        a.cigar = ((0, 10),)
        # Make the T (index 5) have low quality
        quals = [30] * 10
        quals[5] = 5
        a.query_qualities = quals
        outf.write(a)

        # Duplicate (should be ignored if filter enabled)
        a = pysam.AlignedSegment()
        a.query_name = "read_dup"
        a.query_sequence = "AAAAATAAAA"  # Supports SNP T
        a.flag = 1024  # Duplicate
        a.reference_id = 0
        a.reference_start = 95
        a.mapping_quality = 60
        a.cigar = ((0, 10),)
        a.query_qualities = [30] * 10
        outf.write(a)

        # 5. MNP Case: chr1:400, REF=AT, ALT=CG
        # 2 Forward Reads supporting MNP (CG)
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_mnp_fwd_{i}"
            a.query_sequence = "AAAACGAAAA"  # AT -> CG at index 4,5 (pos 400, 401)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 396
            a.mapping_quality = 60
            a.cigar = ((0, 10),)
            a.query_qualities = [30] * 10
            outf.write(a)

        # 1 Reverse Read supporting REF (AT)
        a = pysam.AlignedSegment()
        a.query_name = "read_mnp_ref_rev"
        a.query_sequence = "AAAAATAAAA"
        a.flag = 16
        a.reference_id = 0
        a.reference_start = 396
        a.mapping_quality = 60
        a.cigar = ((0, 10),)
        a.query_qualities = [30] * 10
        outf.write(a)

    # Sort and index the BAM
    sorted_bam = tmp_path / "synthetic.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def test_snp_accuracy(synthetic_bam):
    from gbcms._rs import Variant

    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")

    results = count_bam(
        synthetic_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Expected:
    # Ref Fwd: 5
    # Ref Rev: 3
    # Alt Fwd: 4
    # Alt Rev: 2
    # Ignored: Low MAPQ, Low BaseQ, Duplicate

    assert counts.rd_fwd == 5
    assert counts.rd_rev == 3
    assert counts.ad_fwd == 4
    assert counts.ad_rev == 2
    assert counts.dp == 14  # 5+3+4+2 = 14. (Total reads covering position passing filters)
    # Note: DP might include reads that are neither Ref nor Alt if we had them, but here we only have Ref/Alt.
    # Actually, DP is usually Ref+Alt+Others.
    # Wait, does DP include filtered reads? No, usually raw depth after mapq filter.
    # Let's check implementation: DP increments if mapq passes.
    # But we also check allele status.
    # If !is_ref and !is_alt, we continue. So DP only counts reads that match one of the alleles?
    # Looking at counting.rs:
    # if !is_ref && !is_alt { continue; }
    # counts.dp += 1;
    # So yes, DP only counts reads matching Ref or Alt (or specific alleles if multi-allelic support was better).

    assert counts.rd == 8
    assert counts.ad == 6


def test_insertion_accuracy(synthetic_bam):
    from gbcms._rs import Variant

    variant = Variant(
        chrom="chr1", pos=200, ref_allele="A", alt_allele="AT", variant_type="INSERTION"
    )

    results = count_bam(
        synthetic_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Expected:
    # Alt Fwd: 2
    assert counts.ad_fwd == 2
    assert counts.ad == 2


def test_complex_accuracy(synthetic_bam):
    """
    Test accuracy for complex variants (DelIns, SNP+Indel) using check_complex.
    """
    # bam_path = synthetic_bam
    # 1. Deletion + SNP: REF=AT, ALT=C (Del T, SNP A->C)
    # Read: C (at pos 0), nothing at pos 1.
    # We need to inject a read that has a mismatch and a deletion.
    # This is hard to do with simple pysam.AlignedSegment construction without CIGAR manipulation.
    # But we can construct CIGAR.
    # 2. SNP + Insertion: REF=A, ALT=CT (SNP A->C, Ins T)
    # 3. MNP + Deletion: REF=AAAT, ALT=CG (MNP AA->CG, Del AT)
    # We will create a separate BAM for these complex cases to have full control.
    import pysam

    from gbcms import _rs as gbcms_rs

    complex_bam = "test_complex.bam"

    # Create header
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(complex_bam, "wb", header=header) as outf:
        # Case 1: Del + SNP (REF=AT, ALT=C) at pos 100
        # Read 1: REF (AT) -> Match 2M. Seq=AT
        a = pysam.AlignedSegment()
        a.query_name = "read_del_snp_ref"
        a.query_sequence = "AT"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigar = ((0, 2),)  # 2M
        a.query_qualities = [30, 30]
        outf.write(a)

        # Read 2: ALT (C) -> 1X1D. Seq=C.
        # pysam doesn't support X in CIGAR easily for writing?
        # Actually BAM uses M for match/mismatch usually, but we can use X/=.
        # Let's use M for mismatch (standard BAM).
        # If we use M, we rely on sequence comparison.
        # Ref is AT. Read is C.
        # Alignment: 1M1D.
        # Pos 100: Ref A, Read C -> Mismatch.
        # Pos 101: Ref T, Read Gap -> Deletion.
        b = pysam.AlignedSegment()
        b.query_name = "read_del_snp_alt"
        b.query_sequence = "C"
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 100
        b.mapping_quality = 60
        b.cigar = ((0, 1), (2, 1))  # 1M 1D
        b.query_qualities = [30]
        outf.write(b)

        # Case 2: SNP + Insertion (REF=A, ALT=CT) at pos 200
        # Read 1: REF (A) -> 1M. Seq=A
        c = pysam.AlignedSegment()
        c.query_name = "read_snp_ins_ref"
        c.query_sequence = "A"
        c.flag = 0
        c.reference_id = 0
        c.reference_start = 200
        c.mapping_quality = 60
        c.cigar = ((0, 1),)  # 1M
        c.query_qualities = [30]
        outf.write(c)

        # Read 2: ALT (CT) -> 1M1I. Seq=CT.
        # Ref A, Read C -> Mismatch (1M).
        # Ins T -> 1I.
        d = pysam.AlignedSegment()
        d.query_name = "read_snp_ins_alt"
        d.query_sequence = "CT"
        d.flag = 0
        d.reference_id = 0
        d.reference_start = 200
        d.mapping_quality = 60
        d.cigar = ((0, 1), (1, 1))  # 1M 1I
        d.query_qualities = [30, 30]
        outf.write(d)

    pysam.sort("-o", "test_complex.sorted.bam", complex_bam)
    pysam.index("test_complex.sorted.bam")

    variants = [
        gbcms_rs.Variant("chr1", 100, "AT", "C", "COMPLEX"),
        gbcms_rs.Variant("chr1", 200, "A", "CT", "COMPLEX"),
    ]

    counts = gbcms_rs.count_bam(
        bam_path="test_complex.sorted.bam",
        variants=variants,
        decomposed=[None] * len(variants),
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    # Cleanup
    if os.path.exists(complex_bam):
        os.remove(complex_bam)
    if os.path.exists("test_complex.sorted.bam"):
        os.remove("test_complex.sorted.bam")
    if os.path.exists("test_complex.sorted.bam.bai"):
        os.remove("test_complex.sorted.bam.bai")

    assert len(counts) == 2

    # Case 1: Del + SNP
    # REF read: AT (matches REF AT) -> RD=1
    # ALT read: C (matches ALT C) -> AD=1
    assert counts[0].rd == 1
    assert counts[0].ad == 1

    # Case 2: SNP + Ins
    # REF read: A (matches REF A) -> RD=1
    # ALT read: CT (matches ALT CT) -> AD=1
    assert counts[1].rd == 1
    assert counts[1].ad == 1


def test_deletion_accuracy(synthetic_bam):
    from gbcms._rs import Variant

    variant = Variant(
        chrom="chr1", pos=300, ref_allele="AT", alt_allele="A", variant_type="DELETION"
    )

    results = count_bam(
        synthetic_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Expected:
    # Alt Fwd: 3
    assert counts.ad_fwd == 3
    assert counts.ad == 3


def test_mnp_accuracy(synthetic_bam):
    from gbcms._rs import Variant

    variant = Variant(chrom="chr1", pos=400, ref_allele="AT", alt_allele="CG", variant_type="MNP")

    results = count_bam(
        synthetic_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Expected:
    # Alt Fwd: 2
    # Ref Rev: 1
    assert counts.ad_fwd == 2
    assert counts.rd_rev == 1
    assert counts.ad == 2
    assert counts.rd == 1
