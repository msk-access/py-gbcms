# Docker Guide

Complete guide to using gbcms with Docker.

## Overview

gbcms provides Docker images for:
- **Production use** - Optimized runtime image
- **Testing** - Image with dev dependencies
- **Development** - Local development with docker-compose

## Quick Start

### Pull Pre-built Image

```bash
docker pull mskaccess/gbcms:latest
```

### Run with Docker

```bash
docker run --rm \
    -v /path/to/data:/data \
    mskaccess/gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.txt
```

---

## Production Dockerfile

### Features

✅ **Multi-stage build** - Smaller final image  
✅ **All dependencies** - Includes cyvcf2 for fast VCF parsing  
✅ **Optimized** - Only runtime dependencies in final image  
✅ **Verified** - Installation checked during build  

### What's Included

**System packages**:
- samtools (for BAM/FASTA indexing)
- libhts3 (for cyvcf2)
- All required libraries

**Python packages**:
- gbcms with `[all]` extras
- cyvcf2 (fast VCF parsing)

- Numba (JIT compilation)
- All core dependencies

### Build Locally

```bash
# Build production image
docker build -t gbcms:latest .

# Build with specific tag
docker build -t gbcms:2.0.0 .

# Build with no cache
docker build --no-cache -t gbcms:latest .
```

### Image Size

- **Builder stage**: ~1.5 GB (includes build tools)
- **Final image**: ~800 MB (runtime only)

---

## Usage Examples

### Example 1: Basic VCF Processing

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.txt \
    --thread 4
```

### Example 2: Multiple BAM Files

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --bam sample2:/data/sample2.bam \
    --bam sample3:/data/sample3.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.txt
```

### Example 3: MAF Processing

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam-fof /data/bam_files.txt \
    --maf /data/variants.maf \
    --output /data/counts.maf \
    --omaf
```

### Example 4: With All Options

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run \
    --fasta /data/reference.fa \
    --bam-fof /data/bam_files.txt \
    --vcf /data/variants.vcf.gz \
    --output /data/counts.txt \
    --thread 8 \
    --maq 30 \
    --baq 20 \
    --filter-duplicate \
    --positive-count \
    --fragment-count
```

### Example 5: Validate Files

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    validate files \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf
```

### Example 6: Show Version

```bash
docker run --rm gbcms:latest version
```

### Example 7: Show Help

```bash
docker run --rm gbcms:latest --help
```

---

## Docker Compose

### Using docker-compose.yml

```bash
# Run production container
docker-compose run --rm gbcms count run \
    --fasta /data/reference.fa \
    --bam sample1:/data/sample1.bam \
    --vcf /data/variants.vcf \
    --output /data/counts.txt

# Run tests
docker-compose run --rm test

# Build images
docker-compose build

# Build without cache
docker-compose build --no-cache
```

### Custom docker-compose.yml

```yaml
version: '3.8'

services:
  gbcms:
    image: mskaccess/gbcms:latest
    volumes:
      - ./data:/data
      - ./reference:/reference:ro
    working_dir: /data
    command: >
      count run
      --fasta /reference/genome.fa
      --bam-fof /data/bam_files.txt
      --vcf /data/variants.vcf
      --output /data/counts.txt
      --thread 8
```

---

## Testing Dockerfile

### Purpose

For running tests in a containerized environment.

### Features

✅ **Dev dependencies** - pytest, coverage, linters  
✅ **All features** - cyvcf2, Numba  
✅ **Test fixtures** - Includes test data  

### Build and Run Tests

```bash
# Build test image
docker build -f Dockerfile.test -t gbcms:test .

# Run tests
docker run --rm gbcms:test

# Run specific tests
docker run --rm gbcms:test pytest tests/test_counter.py -v

# Run with coverage report
docker run --rm gbcms:test \
    pytest --cov=gbcms --cov-report=html

# Get coverage report
docker run --rm \
    -v $(pwd)/htmlcov:/app/htmlcov \
    gbcms:test \
    pytest --cov=gbcms --cov-report=html
```

---

## Volume Mounting

### Read-Only Mounts

For reference data that shouldn't be modified:

```bash
docker run --rm \
    -v $(pwd)/reference:/reference:ro \
    -v $(pwd)/bams:/bams:ro \
    -v $(pwd)/output:/output \
    gbcms:latest \
    count run \
    --fasta /reference/genome.fa \
    --bam sample1:/bams/sample1.bam \
    --vcf /bams/variants.vcf \
    --output /output/counts.txt
```

### Multiple Volumes

```bash
docker run --rm \
    -v $(pwd)/reference:/reference:ro \
    -v $(pwd)/bams:/bams:ro \
    -v $(pwd)/variants:/variants:ro \
    -v $(pwd)/output:/output \
    gbcms:latest \
    count run \
    --fasta /reference/genome.fa \
    --bam-fof /bams/bam_files.txt \
    --vcf /variants/variants.vcf \
    --output /output/counts.txt
```

---

## Resource Limits

### CPU Limits

```bash
# Limit to 4 CPUs
docker run --rm \
    --cpus=4 \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run --thread 4 ...

# CPU shares (relative weight)
docker run --rm \
    --cpu-shares=512 \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...
```

### Memory Limits

```bash
# Limit to 8GB RAM
docker run --rm \
    --memory=8g \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...

# Memory with swap
docker run --rm \
    --memory=8g \
    --memory-swap=12g \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...
```

---

## Environment Variables

### Set Environment Variables

```bash
docker run --rm \
    -e NUMBA_NUM_THREADS=4 \
    -e THREADS=16 \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run --thread 16 --backend joblib ...
```

### Using .env File

```bash
# Create .env file
cat > .env << EOF
NUMBA_NUM_THREADS=8
THREADS=16
EOF

# Use with docker run
docker run --rm \
    --env-file .env \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...
```

---

## Troubleshooting

### Issue: Permission Denied

**Problem**: Output files owned by root

**Solution**: Use user mapping

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...
```

### Issue: Out of Memory

**Problem**: Container killed due to OOM

**Solution**: Increase memory limit

```bash
docker run --rm \
    --memory=16g \
    -v $(pwd)/data:/data \
    gbcms:latest \
    count run ...
```

### Issue: Slow Performance

**Problem**: Running slower than expected

**Check**:
1. CPU allocation: `--cpus=N`
2. Thread count: `--thread N`
3. Memory: `--memory=Xg`

### Issue: File Not Found

**Problem**: Cannot find input files

**Solution**: Check volume mounts

```bash
# List files in container
docker run --rm \
    -v $(pwd)/data:/data \
    gbcms:latest \
    bash -c "ls -la /data"
```

### Issue: cyvcf2 Not Working

**Problem**: cyvcf2 import error

**Check**: Image includes cyvcf2

```bash
docker run --rm gbcms:latest \
    python -c "import cyvcf2; print(cyvcf2.__version__)"
```

---

## Best Practices

### 1. Use Specific Tags

```bash
# Good ✅
docker pull mskaccess/gbcms:2.0.0

# Avoid ⚠️
docker pull mskaccess/gbcms:latest
```

### 2. Mount Data as Read-Only When Possible

```bash
docker run --rm \
    -v $(pwd)/reference:/reference:ro \
    -v $(pwd)/bams:/bams:ro \
    -v $(pwd)/output:/output \
    ...
```

### 3. Set Resource Limits

```bash
docker run --rm \
    --cpus=8 \
    --memory=16g \
    ...
```

### 4. Use User Mapping

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    ...
```

### 5. Clean Up

```bash
# Remove stopped containers
docker container prune

# Remove unused images
docker image prune

# Remove all unused data
docker system prune -a
```

---

## CI/CD Integration

### GitHub Actions

```yaml
name: Docker Build

on: [push]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      
      - name: Build Docker image
        run: docker build -t gbcms:latest .
      
      - name: Test Docker image
        run: docker run --rm gbcms:latest version
      
      - name: Run tests
        run: |
          docker build -f Dockerfile.test -t gbcms:test .
          docker run --rm gbcms:test
```

### GitLab CI

```yaml
docker-build:
  stage: build
  script:
    - docker build -t gbcms:latest .
    - docker run --rm gbcms:latest version
  
docker-test:
  stage: test
  script:
    - docker build -f Dockerfile.test -t gbcms:test .
    - docker run --rm gbcms:test
```

---

## Publishing Images

### Docker Hub

```bash
# Tag image
docker tag gbcms:latest mskaccess/gbcms:2.0.0
docker tag gbcms:latest mskaccess/gbcms:latest

# Push to Docker Hub
docker push mskaccess/gbcms:2.0.0
docker push mskaccess/gbcms:latest
```

### GitHub Container Registry

```bash
# Tag for GHCR
docker tag gbcms:latest ghcr.io/msk-access/gbcms:2.0.0

# Login to GHCR
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin

# Push
docker push ghcr.io/msk-access/gbcms:2.0.0
```

---

## Summary

### Dockerfile Features

✅ **Multi-stage build** - Optimized size  
✅ **All dependencies** - cyvcf2, Numba  
✅ **System tools** - samtools included  
✅ **Verified** - Installation checked  
✅ **Documented** - Clear labels  

### Quick Commands

```bash
# Build
docker build -t gbcms:latest .

# Run
docker run --rm -v $(pwd)/data:/data gbcms:latest count run ...

# Test
docker build -f Dockerfile.test -t gbcms:test .
docker run --rm gbcms:test

# Compose
docker-compose run --rm gbcms count run ...
```

### Best Practices

1. Use specific image tags
2. Mount data as read-only when possible
3. Set resource limits
4. Use user mapping for file permissions
5. Clean up unused containers/images

**Docker images are production-ready and include all features!** 🐳
