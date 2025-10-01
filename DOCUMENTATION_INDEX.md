# Documentation Index

Complete index of all GetBaseCounts documentation.

## 📍 GitHub Visibility

**This file (`DOCUMENTATION_INDEX.md`) is in the root directory and will be visible on GitHub.**

However, GitHub only **automatically displays** `README.md` on the repository homepage. To access this index:
- Click on `DOCUMENTATION_INDEX.md` in the file list
- Or link to it from README.md (already done!)

## 📚 Documentation Structure

All documentation is organized in the `docs/` directory following GitBook structure.

### Root Directory Files (Visible on GitHub)
- `README.md` - Main project README (auto-displayed by GitHub)
- `CONTRIBUTING.md` - Contribution guidelines (GitHub links to this)
- `DOCUMENTATION_INDEX.md` - This file (quick reference)
- `LICENSE` - License file

### docs/ Directory (Full Documentation)
All detailed documentation is in `docs/` - see below for complete listing.

### Main Entry Point
- **[docs/README.md](docs/README.md)** - Documentation homepage
- **[docs/SUMMARY.md](docs/SUMMARY.md)** - Table of contents (GitBook)

---

## 🚀 Getting Started (New Users Start Here!)

| Document | Description | Audience |
|----------|-------------|----------|
| [INSTALLATION.md](docs/INSTALLATION.md) | Complete installation guide | Everyone |
| [QUICKSTART.md](docs/QUICKSTART.md) | 5-minute tutorial | New users |
| [CLI_FEATURES.md](docs/CLI_FEATURES.md) | Command-line reference | All users |

---

## 👤 User Guide (Day-to-Day Usage)

| Document | Description | Topics |
|----------|-------------|--------|
| [INPUT_OUTPUT.md](docs/INPUT_OUTPUT.md) | Input/output formats | VCF, MAF, BAM, FASTA, count types |
| [QUALITY_FILTERING.md](docs/QUALITY_FILTERING.md) | Quality filters | Mapping quality, base quality, filters |
| [PERFORMANCE_TUNING.md](docs/PERFORMANCE_TUNING.md) | Optimize performance | Threading, backends, profiling |

---

## 🚀 Advanced Features (Power Users)

| Document | Description | Performance Gain |
|----------|-------------|------------------|
| [ADVANCED_FEATURES.md](docs/ADVANCED_FEATURES.md) | Overview of all features | - |
| [PYDANTIC_GUIDE.md](docs/PYDANTIC_GUIDE.md) | Type safety & validation | Error prevention |
| [CYVCF2_SUPPORT.md](docs/CYVCF2_SUPPORT.md) | Fast VCF parsing | 10-100x faster |
| [NUMBA_GUIDE.md](docs/NUMBA_GUIDE.md) | JIT compilation | 50-100x faster |
| [PARALLELIZATION_GUIDE.md](docs/PARALLELIZATION_GUIDE.md) | Multi-core processing | Linear scaling |
| [RAY_GUIDE.md](docs/RAY_GUIDE.md) | Distributed computing | Cluster scaling |

---

## 🔧 Technical Documentation (Understanding the Code)

| Document | Description | For |
|----------|-------------|-----|
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | Module relationships & data flow | Developers, advanced users |
| [MODULE_GUIDE.md](docs/MODULE_GUIDE.md) | Detailed module descriptions | Developers |
| [COUNTING_ALGORITHMS.md](docs/COUNTING_ALGORITHMS.md) | How counting works | Algorithm developers |
| [GENERIC_COUNTING.md](docs/GENERIC_COUNTING.md) | Generic counting algorithm | Complex variant users |
| [FRAGMENT_COUNTING.md](docs/FRAGMENT_COUNTING.md) | Fragment-level counting | Strand bias analysis |
| [INSERTION_AND_FRAGMENT_ANALYSIS.md](docs/INSERTION_AND_FRAGMENT_ANALYSIS.md) | Detailed analysis | Troubleshooting |

---

## 📖 Reference (Look-Up Information)

| Document | Description | Use When |
|----------|-------------|----------|
| [CPP_FEATURE_COMPARISON.md](docs/CPP_FEATURE_COMPARISON.md) | C++ vs Python comparison | Migrating from C++ |
| [API_REFERENCE.md](docs/API_REFERENCE.md) | Python API documentation | Using as library |
| [CONFIGURATION.md](docs/CONFIGURATION.md) | All config options | Fine-tuning |
| [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | Common issues & solutions | Debugging |
| [FAQ.md](docs/FAQ.md) | Frequently asked questions | Quick answers |
| [GLOSSARY.md](docs/GLOSSARY.md) | Term definitions | Understanding terminology |

## 🐳 Docker & Deployment

| Document | Description | Use When |
|----------|-------------|----------|
| [DOCKER_GUIDE.md](docs/DOCKER_GUIDE.md) | Complete Docker guide | Using Docker |
| [DOCKER_SUMMARY.md](docs/DOCKER_SUMMARY.md) | Docker configuration summary | Quick reference |
| [GITHUB_VISIBILITY.md](docs/GITHUB_VISIBILITY.md) | How docs appear on GitHub | Understanding visibility |
| [DOCUMENTATION_ORGANIZATION.md](docs/DOCUMENTATION_ORGANIZATION.md) | Doc structure | Finding docs |

---

## 💻 Development (Contributors)

| Document | Description | For |
|----------|-------------|-----|
| [CONTRIBUTING.md](CONTRIBUTING.md) | Contribution guidelines | Contributors |
| [DEVELOPMENT.md](docs/DEVELOPMENT.md) | Dev environment setup | Developers |
| [TESTING.md](docs/TESTING.md) | Testing guide | Developers |
| [RELEASE.md](docs/RELEASE.md) | Release process | Maintainers |

---

## 📝 Additional Resources

| Document | Description |
|----------|-------------|
| [FINAL_REVIEW.md](FINAL_REVIEW.md) | Complete project review |
| [PACKAGE_REVIEW.md](PACKAGE_REVIEW.md) | Package status summary |
| [IMPLEMENTATION_SUMMARY.md](IMPLEMENTATION_SUMMARY.md) | Technical implementation details |
| [COMPLETE_FEATURES_SUMMARY.md](COMPLETE_FEATURES_SUMMARY.md) | All features overview |
| [CYVCF2_IMPLEMENTATION_SUMMARY.md](docs/CYVCF2_IMPLEMENTATION_SUMMARY.md) | cyvcf2 implementation |

---

## 🗺️ Documentation Roadmap

### By User Type

#### **New User** (Never used GetBaseCounts)
1. Start: [INSTALLATION.md](docs/INSTALLATION.md)
2. Then: [QUICKSTART.md](docs/QUICKSTART.md)
3. Reference: [CLI_FEATURES.md](docs/CLI_FEATURES.md)
4. Help: [FAQ.md](docs/FAQ.md)

#### **Regular User** (Daily usage)
1. Reference: [INPUT_OUTPUT.md](docs/INPUT_OUTPUT.md)
2. Optimize: [PERFORMANCE_TUNING.md](docs/PERFORMANCE_TUNING.md)
3. Troubleshoot: [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)
4. Advanced: [ADVANCED_FEATURES.md](docs/ADVANCED_FEATURES.md)

#### **Power User** (Large-scale analysis)
1. Speed: [CYVCF2_SUPPORT.md](docs/CYVCF2_SUPPORT.md)
2. Performance: [NUMBA_GUIDE.md](docs/NUMBA_GUIDE.md)
3. Scale: [RAY_GUIDE.md](docs/RAY_GUIDE.md)
4. Optimize: [PARALLELIZATION_GUIDE.md](docs/PARALLELIZATION_GUIDE.md)

#### **Developer** (Contributing code)
1. Setup: [DEVELOPMENT.md](docs/DEVELOPMENT.md)
2. Architecture: [ARCHITECTURE.md](docs/ARCHITECTURE.md)
3. Guidelines: [CONTRIBUTING.md](CONTRIBUTING.md)
4. Testing: [TESTING.md](docs/TESTING.md)

#### **Migrating from C++**
1. Compare: [CPP_FEATURE_COMPARISON.md](docs/CPP_FEATURE_COMPARISON.md)
2. Install: [INSTALLATION.md](docs/INSTALLATION.md)
3. Differences: [FAQ.md](docs/FAQ.md) (see "Comparison with C++")

---

## 📊 Documentation Statistics

| Category | Documents | Total Pages |
|----------|-----------|-------------|
| Getting Started | 3 | ~50 |
| User Guide | 3 | ~75 |
| Advanced Features | 6 | ~150 |
| Technical | 6 | ~200 |
| Reference | 6 | ~100 |
| Development | 4 | ~50 |
| **Total** | **28** | **~625** |

---

## 🔍 Quick Search

### By Topic

**Installation**
- [INSTALLATION.md](docs/INSTALLATION.md)
- [DEVELOPMENT.md](docs/DEVELOPMENT.md)

**Performance**
- [PERFORMANCE_TUNING.md](docs/PERFORMANCE_TUNING.md)
- [CYVCF2_SUPPORT.md](docs/CYVCF2_SUPPORT.md)
- [NUMBA_GUIDE.md](docs/NUMBA_GUIDE.md)
- [PARALLELIZATION_GUIDE.md](docs/PARALLELIZATION_GUIDE.md)

**Counting**
- [COUNTING_ALGORITHMS.md](docs/COUNTING_ALGORITHMS.md)
- [GENERIC_COUNTING.md](docs/GENERIC_COUNTING.md)
- [FRAGMENT_COUNTING.md](docs/FRAGMENT_COUNTING.md)

**File Formats**
- [INPUT_OUTPUT.md](docs/INPUT_OUTPUT.md)
- [CLI_FEATURES.md](docs/CLI_FEATURES.md)

**Troubleshooting**
- [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md)
- [FAQ.md](docs/FAQ.md)

---

## 📱 Documentation Formats

### Online (GitBook)
- Hosted documentation with search and navigation
- Access via: [GitBook URL]

### Markdown (GitHub)
- All docs available in `docs/` directory
- Browse on GitHub with full formatting

### PDF (Generated)
- Complete documentation as single PDF
- Generate with: `make docs-pdf`

---

## 🆘 Getting Help

Can't find what you're looking for?

1. **Search**: Use GitHub search or GitBook search
2. **FAQ**: Check [FAQ.md](docs/FAQ.md)
3. **Issues**: Search [GitHub Issues](https://github.com/msk-access/getbasecounts/issues)
4. **Ask**: Open a [GitHub Discussion](https://github.com/msk-access/getbasecounts/discussions)
5. **Email**: access@mskcc.org

---

## 🔄 Documentation Updates

Documentation is continuously updated. See [CHANGELOG.md](docs/CHANGELOG.md) for recent changes.

**Last Updated**: 2025-10-01

**Version**: 2.0.0

---

**Start exploring**: [docs/README.md](docs/README.md) | [Quick Start](docs/QUICKSTART.md) | [FAQ](docs/FAQ.md)
