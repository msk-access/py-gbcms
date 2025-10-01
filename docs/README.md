# GetBaseCounts Documentation

Welcome to the GetBaseCounts documentation! This guide will help you get started with calculating base counts in BAM files at variant positions.

## What is GetBaseCounts?

GetBaseCounts is a high-performance Python tool for calculating base counts in multiple BAM files at variant positions specified in VCF or MAF files. It's a modern reimplementation of the C++ GetBaseCountsMultiSample tool with enhanced features and performance.

## Key Features

- 🚀 **High Performance**: 50-100x faster with Numba JIT compilation
- ⚡ **Fast VCF Parsing**: 100x faster VCF loading with cyvcf2
- 🔒 **Type Safety**: Runtime validation with Pydantic
- 🔄 **Parallelization**: Multi-core with joblib, distributed with Ray
- 🎨 **Beautiful CLI**: Rich terminal output with progress bars
- 📊 **Multiple Formats**: VCF and MAF input/output
- 🐳 **Docker Support**: Containerized deployment
- 🧪 **Well Tested**: >85% code coverage

## Quick Links

### For New Users
- [Installation Guide](INSTALLATION.md) - Get up and running
- [Quick Start](QUICKSTART.md) - 5-minute tutorial
- [CLI Features](CLI_FEATURES.md) - Command-line reference

### For Advanced Users
- [Advanced Features](ADVANCED_FEATURES.md) - Pydantic, Numba, Ray
- [Performance Tuning](PERFORMANCE_TUNING.md) - Optimize for your workload
- [Architecture](ARCHITECTURE.md) - Understanding the codebase

### For Developers
- [Contributing Guide](../CONTRIBUTING.md) - How to contribute
- [Development Setup](DEVELOPMENT.md) - Set up dev environment
- [API Reference](API_REFERENCE.md) - Python API documentation

## Documentation Structure

This documentation is organized into several sections:

### 📖 Getting Started
Learn the basics and get GetBaseCounts installed and running.

### 👤 User Guide
Detailed guides for using GetBaseCounts in your workflows.

### 🚀 Advanced Features
Deep dives into performance optimization and advanced capabilities.

### 🔧 Technical Documentation
Architecture, algorithms, and implementation details.

### 📚 Reference
Complete reference materials and comparisons.

### 💻 Development
Guides for contributing and developing GetBaseCounts.

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/msk-access/getbasecounts/issues)
- **Discussions**: [GitHub Discussions](https://github.com/msk-access/getbasecounts/discussions)
- **Email**: access@mskcc.org

## Quick Example

```bash
# Install
uv pip install "getbasecounts[all]"

# Run
getbasecounts count run \
    --fasta reference.fa \
    --bam sample1:sample1.bam \
    --vcf variants.vcf \
    --output counts.txt \
    --thread 8
```

## Version

Current version: 2.0.0

## License

Apache License 2.0 - See [LICENSE](../LICENSE) for details.

---

Ready to get started? Head to the [Installation Guide](INSTALLATION.md)!
