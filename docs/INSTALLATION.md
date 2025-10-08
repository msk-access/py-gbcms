# Installation Guide

Complete installation and setup guide for py-gbcms.

## Installation Methods

### Using uv (recommended)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install the package
uv pip install py-gbcms
```

### Using pip

```bash
# Basic installation
pip install py-gbcms

# With fast VCF parsing (10-100x faster)
pip install "py-gbcms[fast]"

# With Ray support for distributed computing
pip install "py-gbcms[ray]"

# With all optional features (fast + Ray)
pip install "py-gbcms[all]"
```

### From source

```bash
git clone https://github.com/msk-access/gbcms.git
cd gbcms

# For production use
uv pip install .

# For development (includes scipy-stubs for type checking)
uv pip install -e ".[dev]"
```

### Using Docker

```bash
docker pull mskaccess/gbcms:latest

# Run the container
docker run --rm \
    -v $(pwd)/data:/data \
    mskaccess/gbcms:latest \
    gbcms count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/output.txt
```

## Verification

### 1. Quick Check

```bash
gbcms version
```

Expected output:
```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                              Version Info                               ┃
┃                                                                          ┃
┃                           py-gbcms                                      ┃
┃                        Version: 2.0.0                                   ┃
┃          Python implementation of gbcms (gbcms)     ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

### 2. Comprehensive Verification

```bash
python3 scripts/verify_installation.py
```

This checks:
- All dependencies installed
- Modules importable
- CLI accessible

### 3. Run Tests

```bash
# Unit tests
pytest

# With coverage
pytest --cov=gbcms

# End-to-end workflow tests
make test-workflows

# Type checking (requires scipy-stubs)
mypy src/
```

## Complete Setup

Run the complete setup and test script:

```bash
make setup
```

This will:
1. Check Python version
2. Install/verify uv
3. Install py-gbcms with all dependencies (including scipy-stubs for type checking)
4. Verify installation
5. Check CLI
6. Run unit tests
7. Run VCF workflow test
8. Run MAF workflow test
