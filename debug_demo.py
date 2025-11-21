import logging
import os

import pysam
from rich.logging import RichHandler

import gbcms_rs

# Setup Logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True)],
)

log = logging.getLogger("rich")
log.info("Starting Debug Demo")

# Create dummy BAM
bam_path = "debug_demo.bam"
header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "read1"
    a.query_sequence = "A" * 100
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 100
    a.mapping_quality = 60
    a.cigar = ((0, 100),)
    outf.write(a)

pysam.index(bam_path)

# Define Variant
# This is a SNP, but let's label it COMPLEX to trigger the check_complex debug logs
variants = [gbcms_rs.Variant("chr1", 150, "A", "T", "COMPLEX")]

# Run Counting
try:
    log.info("Calling Rust engine...")
    gbcms_rs.count_bam(bam_path, variants, 0, 0, True, True, True)
    log.info("Rust engine finished.")
except Exception:
    log.exception("Failed")
finally:
    if os.path.exists(bam_path):
        os.remove(bam_path)
    if os.path.exists(bam_path + ".bai"):
        os.remove(bam_path + ".bai")
