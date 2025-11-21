"""
Test preservation of MAF input columns in output.
"""

import csv

from gbcms.io.input import MafReader
from gbcms.io.output import MafWriter


# Mock counts object
class MockCounts:
    def __init__(self):
        self.dp = 100
        self.rd = 80
        self.ad = 20
        self.dp_fwd = 50
        self.dp_rev = 50
        self.rd_fwd = 40
        self.rd_rev = 40
        self.ad_fwd = 10
        self.ad_rev = 10
        self.dpf = 50
        self.rdf = 40
        self.adf = 10
        self.rdf_fwd = 20
        self.rdf_rev = 20
        self.adf_fwd = 5
        self.adf_rev = 5
        self.sb_pval = 1.0
        self.sb_or = 1.0
        self.fsb_pval = 1.0
        self.fsb_or = 1.0


def test_maf_column_preservation(tmp_path):
    # 1. Create a dummy MAF input with some metadata columns
    input_maf = tmp_path / "input.maf"
    with open(input_maf, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                "Hugo_Symbol",
                "Entrez_Gene_Id",
                "Extra_Column",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Chromosome": "chr1",
                "Start_Position": "100",
                "End_Position": "100",
                "Reference_Allele": "A",
                "Tumor_Seq_Allele2": "T",
                "Hugo_Symbol": "TP53",
                "Entrez_Gene_Id": "7157",
                "Extra_Column": "ShouldBeIgnored",
            }
        )

    # 2. Read the MAF
    reader = MafReader(input_maf)
    variants = list(reader)
    assert len(variants) == 1
    variant = variants[0]

    # Verify metadata was captured
    assert variant.metadata["Hugo_Symbol"] == "TP53"
    assert variant.metadata["Entrez_Gene_Id"] == "7157"
    assert variant.metadata["Extra_Column"] == "ShouldBeIgnored"

    # 3. Write the MAF
    output_maf = tmp_path / "output.maf"
    writer = MafWriter(output_maf)
    writer.write(variant, MockCounts())
    writer.close()

    # 4. Verify output
    with open(output_maf) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

        # Check preserved fields
        assert row["Hugo_Symbol"] == "TP53"
        assert row["Entrez_Gene_Id"] == "7157"

        # Check that calculated fields are present
        assert row["t_total_count"] == "100"

        # Check that extra column is NOT present (since it's not in fieldnames and we used ignore)
        assert "Extra_Column" not in row
