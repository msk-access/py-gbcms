#!/bin/bash
set -e

# Define Paths
DOWNLOADS="$HOME/Downloads"
REF_FASTA="/Users/shahr2/Library/CloudStorage/GoogleDrive-rons.shah@gmail.com/My Drive/Mac_Sep_2024_Documents/test_reference/reference/versions/hg19/Homo_sapiens_assembly19.fasta"

BAM1="$DOWNLOADS/C-6WTKCL-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam"
BAM2="$DOWNLOADS/C-6WTKCL-N001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam"
VCF="$DOWNLOADS/C-U1MJT8-L001-d.DONOR22-TP.combined-variants_anno.vcf"
MAF="$DOWNLOADS/C-000900-L002-d.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots_fillout_filtered.maf"

# Define Output Directories
OUT_VCF_VCF="$DOWNLOADS/gbcms_test_vcf_vcf"
OUT_VCF_MAF="$DOWNLOADS/gbcms_test_vcf_maf"
OUT_MAF_MAF="$DOWNLOADS/gbcms_test_maf_maf"
OUT_MAF_VCF="$DOWNLOADS/gbcms_test_maf_vcf"
OUT_MULTI="$DOWNLOADS/gbcms_test_multi"

# Ensure Output Directories Exist (gbcms creates them, but good to be clean)
rm -rf "$OUT_VCF_VCF" "$OUT_VCF_MAF" "$OUT_MAF_MAF" "$OUT_MAF_VCF" "$OUT_MULTI"

echo "========================================================"
echo "Starting Real Data Verification"
echo "Reference: $REF_FASTA"
echo "BAM1: $BAM1"
echo "VCF: $VCF"
echo "MAF: $MAF"
echo "========================================================"

# 1. VCF Input -> VCF Output
echo ""
echo "[TEST 1] VCF Input -> VCF Output"
python -m gbcms.cli run \
    --variants "$VCF" \
    --bam "$BAM1" \
    --fasta "$REF_FASTA" \
    --output-dir "$OUT_VCF_VCF" \
    --format vcf \
    --verbose

# 2. VCF Input -> MAF Output
echo ""
echo "[TEST 2] VCF Input -> MAF Output"
python -m gbcms.cli run \
    --variants "$VCF" \
    --bam "$BAM1" \
    --fasta "$REF_FASTA" \
    --output-dir "$OUT_VCF_MAF" \
    --format maf \
    --verbose

# 3. MAF Input -> MAF Output
echo ""
echo "[TEST 3] MAF Input -> MAF Output"
python -m gbcms.cli run \
    --variants "$MAF" \
    --bam "$BAM1" \
    --fasta "$REF_FASTA" \
    --output-dir "$OUT_MAF_MAF" \
    --format maf \
    --verbose

# 4. MAF Input -> VCF Output
echo ""
echo "[TEST 4] MAF Input -> VCF Output"
python -m gbcms.cli run \
    --variants "$MAF" \
    --bam "$BAM1" \
    --fasta "$REF_FASTA" \
    --output-dir "$OUT_MAF_VCF" \
    --format vcf \
    --verbose

# 5. Multi-BAM Test (VCF -> MAF)
echo ""
echo "[TEST 5] Multi-BAM Test (VCF -> MAF)"
python -m gbcms.cli run \
    --variants "$VCF" \
    --bam "$BAM1" \
    --bam "$BAM2" \
    --fasta "$REF_FASTA" \
    --output-dir "$OUT_MULTI" \
    --format maf \
    --verbose

echo ""
echo "========================================================"
echo "All tests completed successfully!"
echo "Outputs are in:"
echo "  - $OUT_VCF_VCF"
echo "  - $OUT_VCF_MAF"
echo "  - $OUT_MAF_MAF"
echo "  - $OUT_MAF_VCF"
echo "  - $OUT_MULTI"
echo "========================================================"
