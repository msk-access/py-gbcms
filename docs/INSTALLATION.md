# Installation Guide

This guide covers all methods for installing py-gbcms.

## Quick Install (Recommended)

### PyPI Installation

Install the latest stable release:

```bash
pip install py-gbcms
```

Verify installation:

```bash
gbcms --version
```

### Docker

Pull the official Docker image:

```bash
docker pull ghcr.io/msk-access/py-gbcms:2.1.0
```

Run via Docker:

```bash
docker run --rm ghcr.io/msk-access/py-gbcms:2.1.0 gbcms --help
```

---

## Requirements

### System Requirements

- **Python**: 3.10 or higher
- **Operating Systems**: Linux, macOS, Windows (with WSL2)
- **Memory**: Minimum 4GB RAM (8GB+ recommended for large BAM files)
- **Disk**: Varies by BAM file size

### For Source Installation

- **Rust toolchain**: 1.70 or higher (for building from source)
- **Build tools**: `gcc`, `make` (Linux/macOS)

### For Nextflow Workflow

- **Nextflow**: 21.10.3 or higher
- **Container runtime**: Docker or Singularity
- **HPC scheduler**: SLURM (optional, for cluster execution)

---

## Installation Methods

### 1. PyPI (Stable Release)

**Best for:** Most users, production use

```bash
pip install py-gbcms
```

**Benefits:**
- Pre-compiled wheels for major platforms
- No Rust toolchain required
- Fastest installation

### 2. From Source (Development)

**Best for:** Developers, contributing code, bleeding-edge features

#### Prerequisites

Install Rust toolchain:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

#### Installation Steps

```bash
# Clone repository
git clone https://github.com/msk-access/py-gbcms.git
cd py-gbcms

# Install with uv (recommended)
curl -LsSf https://astral.sh/uv/install.sh | sh
uv pip install -e ".[dev]"

# OR install with pip
pip install -e ".[dev]"
```

Verify:

```bash
gbcms --version
```

### 3. Docker

**Best for:** Reproducibility, containerized workflows, no local dependencies

#### Pull Image

```bash
# Latest stable release
docker pull ghcr.io/msk-access/py-gbcms:2.1.0

# Or use 'latest' tag
docker pull ghcr.io/msk-access/py-gbcms:latest
```

#### Usage Example

```bash
docker run --rm \
  -v $(pwd):/data \
  ghcr.io/msk-access/py-gbcms:2.1.0 \
  gbcms run \
    --variants /data/variants.vcf \
    --bam /data/sample.bam \
    --fasta /data/reference.fa \
    --output-dir /data/results/
```

**Note:** Use `-v` to mount your data directory into the container.

### 4. Singularity (HPC)

**Best for:** HPC clusters without Docker

```bash
# Pull from Docker registry
singularity pull docker://ghcr.io/msk-access/py-gbcms:2.1.0

# Run
singularity exec py-gbcms_2.1.0.sif gbcms --help
```

---

## Nextflow Workflow Setup

For processing multiple samples in parallel:

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or any directory in your PATH
```

Verify:

```bash
nextflow -version
```

### 2. Install Container Runtime

**Docker (local):**

Follow instructions at https://docs.docker.com/engine/installation/

**Singularity (HPC):**

Singularity is usually pre-installed on HPC systems. Check with:

```bash
singularity --version
```

If not available, contact your system administrator.

### 3. Clone Repository (for Nextflow workflow)

```bash
git clone https://github.com/msk-access/py-gbcms.git
cd py-gbcms
```

The Nextflow workflow is in the `nextflow/` directory.

---

## Verification

### Test Installation

Run a simple command:

```bash
gbcms --help
```

Expected output:

```
Usage: gbcms [OPTIONS] COMMAND [ARGS]...

  py-gbcms: Get Base Counts Multi-Sample

Options:
  --version  Show version and exit
  --help     Show this message and exit

Commands:
  run  Run variant counting
```

### Check Version

```bash
gbcms --version
```

Should show: `2.1.0` (or current version)

---

## Troubleshooting

### Python Version Issues

**Error:** `Python 3.10+ required`

**Solution:**

```bash
# Check Python version
python --version

# Install Python 3.10+ (Ubuntu/Debian)
sudo apt update
sudo apt install python3.10 python3.10-pip

# Use pyenv (any OS)
pyenv install 3.10.0
pyenv global 3.10.0
```

### Rust Toolchain Not Found

**Error:** `cargo not found` or `rustc not found`

**Solution:**

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Reload shell
source $HOME/.cargo/env

# Verify
rustc --version
```

### Docker Permission Denied

**Error:** `permission denied while trying to connect to the Docker daemon`

**Solution:**

```bash
# Add user to docker group (Linux)
sudo usermod -aG docker $USER
newgrp docker

# Or use sudo
sudo docker pull ghcr.io/msk-access/py-gbcms:2.1.0
```

### Module Import Errors

**Error:** `ModuleNotFoundError: No module named 'gbcms'`

**Solution:**

```bash
# Reinstall
pip uninstall py-gbcms
pip install py-gbcms

# Or force reinstall
pip install --force-reinstall py-gbcms
```

### BAM Index Not Found

**Error:** `BAI index not found`

**Solution:**

```bash
# Create BAM index
samtools index sample.bam
```

This creates `sample.bam.bai` in the same directory.

---

## Upgrading

### Upgrade from PyPI

```bash
pip install --upgrade py-gbcms
```

### Upgrade Docker Image

```bash
docker pull ghcr.io/msk-access/py-gbcms:latest
```

### Upgrade from Source

```bash
cd py-gbcms
git pull origin main
pip install -e ".[dev]"
```

---

## Uninstallation

### PyPI

```bash
pip uninstall py-gbcms
```

### Docker

```bash
docker rmi ghcr.io/msk-access/py-gbcms:2.1.0
```

### Singularity

```bash
rm py-gbcms_2.1.0.sif
```

---

## Next Steps

After installation:

- **CLI Users:** See [CLI Quick Start](quick-start.md)
- **Nextflow Users:** See [Nextflow Workflow Guide](NEXTFLOW.md)
- **Choosing a method:** See [Usage Patterns](WORKFLOWS.md)
