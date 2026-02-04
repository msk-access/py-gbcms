# Installation

## Quick Install

=== "PyPI (Recommended)"
    ```bash
    pip install py-gbcms
    gbcms --version
    ```

=== "Docker"
    ```bash
    docker pull ghcr.io/msk-access/py-gbcms:2.2.0
    docker run --rm ghcr.io/msk-access/py-gbcms:2.2.0 gbcms --help
    ```

=== "From Source"
    ```bash
    git clone https://github.com/msk-access/py-gbcms.git
    cd py-gbcms
    pip install -e ".[dev]"
    ```

---

## Requirements

| Component | Requirement |
|:----------|:------------|
| Python | 3.10+ |
| OS | Linux, macOS, Windows (WSL2) |
| Memory | 4GB+ (8GB for large BAMs) |

### For Nextflow Workflow

- Nextflow 21.10.3+
- Docker or Singularity

---

## Verification

```bash
# Check installation
gbcms --version
# Expected: 2.2.0

# Test help
gbcms --help
```

---

## Docker Usage

```bash
docker run --rm \
  -v $(pwd):/data \
  ghcr.io/msk-access/py-gbcms:2.2.0 \
  gbcms run \
    --variants /data/variants.vcf \
    --bam /data/sample.bam \
    --fasta /data/reference.fa \
    --output-dir /data/results/
```

> **Tip:** Use `-v` to mount your data directory.

---

## Singularity (HPC)

```bash
singularity pull docker://ghcr.io/msk-access/py-gbcms:2.2.0
singularity exec py-gbcms_2.2.0.sif gbcms --help
```

---

## Troubleshooting

### Module Not Found
```bash
pip uninstall py-gbcms && pip install py-gbcms
```

### BAM Index Missing
```bash
samtools index sample.bam
```

### Docker Permission Denied
```bash
sudo usermod -aG docker $USER && newgrp docker
```

---

## Upgrade

```bash
pip install --upgrade py-gbcms
```

---

## Next Steps

- **[CLI Quick Start](quick-start.md)** — Command examples
- **[Nextflow Guide](NEXTFLOW.md)** — HPC pipeline
