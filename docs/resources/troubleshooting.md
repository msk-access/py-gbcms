# Troubleshooting

Common issues and solutions for py-gbcms.

## Installation Issues

### Rust Compilation Fails

**Error:** `error: linker 'cc' not found`

**Solution:** Install build tools:
```bash
# macOS
xcode-select --install

# Ubuntu/Debian
apt-get install build-essential

# CentOS/RHEL
yum groupinstall "Development Tools"
```

### htslib Linking Error

**Error:** `cannot find -lhts`

**Solution:** The Rust backend compiles htslib statically. Ensure you have:
```bash
# Dependencies
apt-get install zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
```

## Runtime Issues

### BAI Not Found

**Error:** `BAI index not found for sample.bam`

**Solution:** Create or verify index:
```bash
samtools index sample.bam
```

### Chromosome Mismatch

**Error:** `Chromosome X not found in reference`

**Solution:** Ensure chromosome naming matches between VCF/MAF, BAM, and FASTA (e.g., `chr1` vs `1`).

## Nextflow Issues

### Container Not Found

**Error:** `Unable to find image 'ghcr.io/msk-access/py-gbcms:latest'`

**Solution:** Pull manually:
```bash
docker pull ghcr.io/msk-access/py-gbcms:latest
# or
singularity pull docker://ghcr.io/msk-access/py-gbcms:latest
```

### Task Metrics Error

**Error:** `Command 'ps' required by nextflow to collect task metrics`

This occurs when running containers that lack `procps`. Use a newer image with procps installed.

## Related

- [Installation](../getting-started/installation.md) — Setup guide
- [Nextflow Parameters](../nextflow/parameters.md) — Resource options
