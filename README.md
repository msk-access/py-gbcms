# py-gbcms

**Complete orientation-aware counting system for genomic variants**

[![Tests](https://github.com/msk-access/py-gbcms/workflows/Tests/badge.svg)](https://github.com/msk-access/py-gbcms/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Features

- 🚀 **High Performance**: Rust-powered core engine with multi-threading
- 🧬 **Complete Variant Support**: SNP, MNP, insertion, deletion, and complex variants (DelIns, SNP+Indel)
- 📊 **Orientation-Aware**: Forward and reverse strand analysis with fragment counting
- 🔬 **Statistical Analysis**: Fisher's exact test for strand bias
- 📁 **Flexible I/O**: VCF and MAF input/output formats
- 🎯 **Quality Filters**: 7 configurable read filtering options

## Installation

**Quick install:**
```bash
pip install py-gbcms
```

**From source (requires Rust):**
```bash
git clone https://github.com/msk-access/py-gbcms.git
cd py-gbcms
pip install .
```

**Docker:**
```bash
docker pull ghcr.io/msk-access/py-gbcms:2.1.0
```

📖 **Full documentation:** https://msk-access.github.io/py-gbcms/

---

## Usage

`py-gbcms` can be used in two ways:

### 🔧 Option 1: Standalone CLI (1-10 samples)

**Best for:** Quick analysis, local processing, direct control

```bash
gbcms run \
    --variants variants.vcf \
    --bam sample1.bam \
    --fasta reference.fa \
    --output-dir results/
```

**Output:** `results/sample1.vcf`

**Learn more:**
- 📘 [CLI Quick Start](https://cmo-ci.gitbook.io/py-gbcms/quick-start)
- 📖 [CLI Reference](https://cmo-ci.gitbook.io/py-gbcms/cli_features)

---

### 🔄 Option 2: Nextflow Workflow (10+ samples, HPC)

**Best for:** Many samples, HPC clusters (SLURM), reproducible pipelines

```bash
nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta reference.fa \
    -profile slurm
```

**Features:**
- ✅ Automatic parallelization across samples
- ✅ SLURM/HPC integration
- ✅ Container support (Docker/Singularity)
- ✅ Resume failed runs

**Learn more:**
- 🔄 [Nextflow Workflow Guide](https://cmo-ci.gitbook.io/py-gbcms/nextflow)
- 📋 [Usage Patterns Comparison](https://cmo-ci.gitbook.io/py-gbcms/workflows)

---

## Which Should I Use?

| Scenario | Recommendation |
|----------|----------------|
| 1-10 samples, local machine | **CLI** |
| 10+ samples, HPC cluster | **Nextflow** |
| Quick ad-hoc analysis | **CLI** |
| Production pipeline | **Nextflow** |
| Need auto-parallelization | **Nextflow** |
| Full manual control | **CLI** |

---

## Quick Examples

### CLI: Single Sample
```bash
gbcms run \
    --variants variants.vcf \
    --bam tumor.bam \
    --fasta hg19.fa \
    --output-dir results/ \
    --threads 4
```

### CLI: Multiple Samples (Sequential)
```bash
gbcms run \
    --variants variants.vcf \
    --bam-list samples.txt \
    --fasta hg19.fa \
    --output-dir results/
```

### Nextflow: Many Samples (Parallel)
```bash
# samplesheet.csv:
# sample,bam,bai
# tumor1,/path/to/tumor1.bam,
# tumor2,/path/to/tumor2.bam,

nextflow run nextflow/main.nf \
    --input samplesheet.csv \
    --variants variants.vcf \
    --fasta hg19.fa \
    --outdir results \
    -profile slurm
```

---

## Documentation

📚 **Full Documentation:** https://cmo-ci.gitbook.io/py-gbcms/

**Quick Links:**
- [Installation](https://cmo-ci.gitbook.io/py-gbcms/installation)
- [CLI Quick Start](https://cmo-ci.gitbook.io/py-gbcms/quick-start)
- [Nextflow Workflow](https://cmo-ci.gitbook.io/py-gbcms/nextflow)
- [CLI Reference](https://cmo-ci.gitbook.io/py-gbcms/cli_features)
- [Input & Output Formats](https://cmo-ci.gitbook.io/py-gbcms/input_output)
- [Architecture](https://cmo-ci.gitbook.io/py-gbcms/architecture)

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines.

To contribute to documentation, see the [`gh-pages` branch](https://github.com/msk-access/py-gbcms/tree/gh-pages).

---

## Citation

If you use `py-gbcms` in your research, please cite:

```
[Citation to be added]
```

---

## License

AGPL-3.0 - see [LICENSE](LICENSE) for details.

---

## Support

- 🐛 **Issues:** https://github.com/msk-access/py-gbcms/issues
- 💬 **Discussions:** https://github.com/msk-access/py-gbcms/discussions
