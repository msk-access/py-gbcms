#!/bin/bash
# Validation script to compare Python implementation against C++ version

set -e

echo "=========================================="
echo "Python vs C++ Validation"
echo "=========================================="
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if C++ binary exists
CPP_BINARY="/Users/shahr2/Documents/Github/GetBaseCountsMultiSample/GetBaseCountsMultiSample"
if [ ! -f "$CPP_BINARY" ]; then
    echo -e "${YELLOW}⚠️  C++ binary not found at: $CPP_BINARY${NC}"
    echo "Please compile the C++ version first:"
    echo "  cd /Users/shahr2/Documents/Github/GetBaseCountsMultiSample"
    echo "  make"
    exit 1
fi

# Create test directory
TEST_DIR=$(mktemp -d)
echo "Test directory: $TEST_DIR"
cd "$TEST_DIR"

# Create test data
echo "Creating test data..."

# Reference FASTA
cat > reference.fa << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
EOF

samtools faidx reference.fa

# VCF file
cat > variants.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=180>
##contig=<ID=chr2,length=180>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	T	.	PASS	.
chr1	20	.	C	G	.	PASS	.
chr1	30	.	AT	A	.	PASS	.
chr1	40	.	G	GC	.	PASS	.
chr2	15	.	G	C	.	PASS	.
EOF

# BAM file
samtools view -bS > sample1.bam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:180
@SQ	SN:chr2	LN:180
read1	99	chr1	5	60	20M	=	50	65	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read2	99	chr1	15	60	20M	=	60	65	CGATCGATCGATCGATCGAT	IIIIIIIIIIIIIIIIIIII
read3	99	chr1	25	60	20M	=	70	65	GATCGATCGATCGATCGATC	IIIIIIIIIIIIIIIIIIII
read4	99	chr1	35	60	20M	=	80	65	TCGATCGATCGATCGATCGA	IIIIIIIIIIIIIIIIIIII
read5	99	chr2	10	60	20M	=	55	65	GCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIII
EOF

samtools index sample1.bam

echo

# Test 1: Basic VCF workflow
echo "=========================================="
echo "Test 1: Basic VCF Workflow"
echo "=========================================="

echo "Running C++ version..."
$CPP_BINARY \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output cpp_output.txt \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo "Running Python version..."
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output py_output.txt \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo
echo "Comparing outputs..."
if diff -q cpp_output.txt py_output.txt > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Test 1 PASSED: Outputs match!${NC}"
    TEST1_PASS=1
else
    echo -e "${YELLOW}⚠️  Test 1: Outputs differ${NC}"
    echo "Showing first 10 lines of each:"
    echo "--- C++ ---"
    head -10 cpp_output.txt
    echo "--- Python ---"
    head -10 py_output.txt
    echo "--- Diff ---"
    diff cpp_output.txt py_output.txt | head -20 || true
    TEST1_PASS=0
fi

echo

# Test 2: With quality filters
echo "=========================================="
echo "Test 2: With Quality Filters"
echo "=========================================="

echo "Running C++ version..."
$CPP_BINARY \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output cpp_filtered.txt \
    --maq 30 \
    --baq 20 \
    --filter_duplicate 1 \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo "Running Python version..."
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output py_filtered.txt \
    --maq 30 \
    --baq 20 \
    --filter-duplicate \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo
echo "Comparing outputs..."
if diff -q cpp_filtered.txt py_filtered.txt > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Test 2 PASSED: Filtered outputs match!${NC}"
    TEST2_PASS=1
else
    echo -e "${YELLOW}⚠️  Test 2: Filtered outputs differ${NC}"
    TEST2_PASS=0
fi

echo

# Test 3: Generic counting
echo "=========================================="
echo "Test 3: Generic Counting"
echo "=========================================="

echo "Running C++ version with --generic_counting..."
$CPP_BINARY \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output cpp_generic.txt \
    --generic_counting \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo "Running Python version with --generic-counting..."
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output py_generic.txt \
    --generic-counting \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo
echo "Comparing outputs..."
if diff -q cpp_generic.txt py_generic.txt > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Test 3 PASSED: Generic counting outputs match!${NC}"
    TEST3_PASS=1
else
    echo -e "${YELLOW}⚠️  Test 3: Generic counting outputs differ${NC}"
    TEST3_PASS=0
fi

echo

# Test 4: Fragment counting
echo "=========================================="
echo "Test 4: Fragment Counting"
echo "=========================================="

echo "Running C++ version with --fragment_count..."
$CPP_BINARY \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output cpp_fragment.txt \
    --fragment_count 1 \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo "Running Python version with --fragment-count..."
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output py_fragment.txt \
    --fragment-count \
    --thread 1 \
    2>&1 | grep -v "Warning" || true

echo
echo "Comparing outputs..."
if diff -q cpp_fragment.txt py_fragment.txt > /dev/null 2>&1; then
    echo -e "${GREEN}✅ Test 4 PASSED: Fragment counting outputs match!${NC}"
    TEST4_PASS=1
else
    echo -e "${YELLOW}⚠️  Test 4: Fragment counting outputs differ${NC}"
    TEST4_PASS=0
fi

echo

# Summary
echo "=========================================="
echo "Summary"
echo "=========================================="

TOTAL_TESTS=4
PASSED_TESTS=$((TEST1_PASS + TEST2_PASS + TEST3_PASS + TEST4_PASS))

echo "Tests passed: $PASSED_TESTS / $TOTAL_TESTS"
echo

if [ $PASSED_TESTS -eq $TOTAL_TESTS ]; then
    echo -e "${GREEN}✅ All tests PASSED!${NC}"
    echo "Python implementation matches C++ implementation."
    EXIT_CODE=0
else
    echo -e "${YELLOW}⚠️  Some tests failed or showed differences${NC}"
    echo "This may be due to:"
    echo "  - Floating point precision differences"
    echo "  - Different random tie-breaking"
    echo "  - Minor implementation differences"
    echo
    echo "Review the differences above to determine if they are acceptable."
    EXIT_CODE=1
fi

# Cleanup
echo
echo "Cleaning up..."
cd /
rm -rf "$TEST_DIR"

echo
echo "=========================================="
echo "Validation complete"
echo "=========================================="

exit $EXIT_CODE
