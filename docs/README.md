# 📚 gbcms Documentation Hub

Welcome to the comprehensive documentation for **gbcms** - the Python implementation of gbcms.

## 🚀 Quick Start

| Action | Command | Documentation |
|--------|---------|---------------|
| **Install** | `uv pip install "gbcms[all]"` | [Installation Guide](INSTALLATION.md) |
| **Basic Usage** | `gbcms count run --fasta ref.fa --bam sample1.bam --vcf variants.vcf --output counts.txt` | [Quick Start](QUICKSTART.md) |
| **Docker** | `docker run ghcr.io/msk-access/gbcms:latest gbcms count run --omaf --fasta /data/ref.fa --bam sample1:/data/sample1.bam --vcf /data/variants.vcf --output /data/counts.maf` | [Docker Guide](DOCKER_GUIDE.md) |

## 📖 Documentation Sections

### 🔰 **Getting Started**
- **[Installation Guide](INSTALLATION.md)** - Complete setup instructions for different environments
- **[Quick Start Tutorial](QUICKSTART.md)** - 5-minute guide to process your first samples
- **[CLI Features](CLI_FEATURES.md)** - Complete command-line interface reference

### 👥 **User Guides**
- **[Input & Output Formats](INPUT_OUTPUT.md)** - Supported file formats and specifications

### ⚡ **Advanced Features**
- **[Advanced Features](ADVANCED_FEATURES.md)** - Pydantic, Numba, Ray deep dive
- **[Parallelization Guide](PARALLELIZATION_GUIDE.md)** - When and how to use joblib vs Ray
- **[Fast VCF Parsing (cyvcf2)](CYVCF2_SUPPORT.md)** - High-performance VCF processing

### 🏗️ **Architecture & Development**
- **[Architecture Overview](ARCHITECTURE.md)** - Module relationships and design decisions
- **[C++ Feature Comparison](CPP_FEATURE_COMPARISON.md)** - How Python version compares to C++
- **[FAQ](FAQ.md)** - Frequently asked questions
- **[Troubleshooting](TROUBLESHOOTING.md)** - Common issues and solutions

### 🐳 **Docker & Deployment**
- **[Docker Guide](DOCKER_GUIDE.md)** - Complete containerization guide
- **[Docker Summary](DOCKER_SUMMARY.md)** - Quick reference for Docker usage

## 🛠️ Development

### **For Contributors**
- **[CONTRIBUTING.md](../../../CONTRIBUTING.md)** - Development workflow and guidelines
- **[Package Structure](PACKAGE_STRUCTURE.md)** - Understanding the codebase organization
- **[Testing Guide](TESTING_GUIDE.md)** - Testing strategy and coverage

### **Code Quality**
```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/

# Type checking
mypy src/

# Run tests with coverage
pytest --cov=gbcms --cov-report=html
```

## 📊 Performance Benchmarks

| Feature | C++ Version | Python (Basic) | Python (Optimized) |
|---------|-------------|----------------|-------------------|
| Speed | ~1x | ~0.8-1.2x | **~2-5x** |
| Memory | Baseline | ~1.2x | ~1.5x |
| Multi-threading | OpenMP | concurrent.futures | **joblib/Ray** |
| Scalability | Single machine | Single machine | **Multi-node clusters** |

**Python 3.11+ shows significant improvements with Numba JIT compilation.**

## 🔗 Key Links

- **🐛 [Report Issues](https://github.com/msk-access/py-gbcms/issues)**
- **💬 [Discussions](https://github.com/msk-access/py-gbcms/discussions)**
- **📧 [Contact](mailto:shahr2@mskcc.org)** - MSK-ACCESS Team
- **📦 [PyPI Package](https://pypi.org/project/py-gbcms/)**
- **🐳 [Docker Hub](https://ghcr.io/msk-access/py-gbcms)**

## 📝 Citation

If you use gbcms in your research, please cite:

```bibtex
@software{gbcms,
  title={gbcms: Python implementation of gbcms},
  author={MSK-ACCESS Team},
  url={https://github.com/msk-access/py-gbcms},
  version={2.0.0}
}
```

## 📄 License

**AGPL-3.0 License** - See [LICENSE](../../../LICENSE) for details.

---

**Built with ❤️ by the MSK-ACCESS Team**

*This is a Python reimplementation of the original C++ tool for enhanced maintainability and extensibility while maintaining high performance through modern Python optimization techniques.*
