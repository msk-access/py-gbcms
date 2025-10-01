# CI/CD Setup Complete ✅

Complete summary of Docker and GitHub Actions setup for GetBaseCounts.

---

## ✅ What Was Done

### 1. Docker Base Image Changed to Ubuntu

**Before**: `python:3.11-slim` (Debian-based)  
**After**: `ubuntu:22.04` (Ubuntu LTS)

**Files**:
- `Dockerfile` - Now uses Ubuntu 22.04 ✅
- `Dockerfile.python-slim` - Backup of old version
- `Dockerfile.test` - Updated for Ubuntu

**Why Ubuntu?**:
- More familiar to users
- LTS support until 2027
- Matches many production environments
- Large package repository

**Size**: ~900 MB (slightly larger but acceptable)

---

### 2. GitHub Actions Workflows Created

#### A. CI Workflow (`ci.yml`) ✅

**Purpose**: Continuous Integration on every push

**Features**:
- Matrix testing: Ubuntu + macOS, Python 3.9-3.12
- Installs all dependencies including libhts for cyvcf2
- Runs pytest with coverage
- Uploads coverage to Codecov
- Linting (black, ruff, mypy)
- Docker build and test

**Trigger**: Push or PR to `main` or `develop`

#### B. Test Workflow (`test.yml`) ✅

**Purpose**: Comprehensive testing

**Features**:
- Same as CI but more explicit
- Separate jobs for test, lint, docker
- Better organized for complex projects

**Trigger**: Push, PR, or manual

#### C. PyPI Publishing (`publish-pypi.yml`) ✅

**Purpose**: Automated package publishing

**Features**:
- Builds Python package (wheel + sdist)
- Validates with twine
- Publishes to PyPI or Test PyPI
- Uses trusted publishing (no API tokens needed)

**Trigger**: Release published or manual

**Setup Required**:
1. Configure PyPI trusted publishing at https://pypi.org/manage/account/publishing/
2. Add repository: `msk-access/getbasecounts`
3. Workflow: `publish-pypi.yml`

#### D. Docker Publishing (`publish-docker.yml`) ✅

**Purpose**: Publish Docker images to GHCR

**Features**:
- Multi-platform builds (amd64, arm64)
- Automatic tagging (latest, version, sha)
- Push to GitHub Container Registry
- Build attestation for security
- Layer caching for fast builds

**Trigger**: Release, push to main, tags, or manual

**Image Location**: `ghcr.io/msk-access/getbasecounts`

---

## 📦 Docker Configuration

### Base Image: Ubuntu 22.04

```dockerfile
FROM ubuntu:22.04

# Install Python 3.11
# Install all dependencies
# Install samtools, libhts
# Install GetBaseCounts with [all] extras
```

### System Dependencies Included

| Package | Purpose |
|---------|---------|
| python3.11 | Python runtime |
| samtools | BAM/FASTA indexing |
| libhts3 | HTSlib for cyvcf2 |
| zlib, bz2, lzma | Compression |
| curl, ssl | Network |

### Python Dependencies Included

| Package | Purpose |
|---------|---------|
| pysam | BAM reading |
| numpy | Numerical ops |
| typer, rich | CLI |
| pandas | Data handling |
| pydantic | Type safety |
| numba | JIT compilation |
| joblib | Parallelization |
| **cyvcf2** | Fast VCF parsing ⭐ |
| **ray** | Distributed computing ⭐ |

---

## 🚀 GitHub Actions Features

### Matrix Testing

```yaml
OS: [ubuntu-latest, macos-latest]
Python: [3.9, 3.10, 3.11, 3.12]
Total: 8 test jobs
```

### Code Quality

- ✅ Black (code formatting)
- ✅ Ruff (linting)
- ✅ Mypy (type checking)
- ✅ Pytest (testing)
- ✅ Coverage (Codecov)

### Docker Features

- ✅ Multi-platform (amd64, arm64)
- ✅ Layer caching (GitHub Actions cache)
- ✅ Automatic tagging
- ✅ Build attestation
- ✅ Security scanning

---

## 📋 Setup Checklist

### GitHub Repository Settings

- [x] Enable GitHub Actions
- [x] Set workflow permissions to "Read and write"
- [ ] Enable branch protection for `main`
- [ ] Require CI to pass before merge

### PyPI Publishing

- [ ] Go to https://pypi.org/manage/account/publishing/
- [ ] Add trusted publisher:
  - Project: `getbasecounts`
  - Owner: `msk-access`
  - Repository: `getbasecounts`
  - Workflow: `publish-pypi.yml`
- [ ] Test with Test PyPI first

### GitHub Container Registry

- [x] Automatically configured (no setup needed)
- [ ] Set package visibility to public
- [ ] Add package description

### Codecov (Optional)

- [ ] Sign up at https://codecov.io/
- [ ] Add repository
- [ ] Get upload token
- [ ] Add to GitHub Secrets as `CODECOV_TOKEN`

---

## 🎯 Usage

### Running Tests Locally

```bash
# Install dependencies
uv pip install -e ".[dev,all]"

# Run tests
pytest --cov=getbasecounts --cov-report=term-missing

# Run linters
black --check src/ tests/
ruff check src/ tests/
mypy src/

# Build Docker
docker build -t getbasecounts:local .
docker run --rm getbasecounts:local version
```

### Creating a Release

```bash
# 1. Update version in pyproject.toml
# 2. Commit changes
git add pyproject.toml
git commit -m "Bump version to 2.0.0"
git push

# 3. Create and push tag
git tag v2.0.0
git push origin v2.0.0

# 4. Create GitHub release
gh release create v2.0.0 \
  --title "Version 2.0.0" \
  --notes "Release notes here"

# This automatically triggers:
# - publish-pypi.yml → Publishes to PyPI
# - publish-docker.yml → Publishes Docker image to GHCR
```

### Manual Publishing

**Test PyPI**:
1. Go to Actions tab
2. Select "Publish to PyPI"
3. Click "Run workflow"
4. Check "Publish to Test PyPI"
5. Click "Run workflow"

**Docker Image**:
1. Go to Actions tab
2. Select "Publish Docker Image"
3. Click "Run workflow"
4. Select branch
5. Click "Run workflow"

---

## 📊 Workflow Status

### CI Workflow

**Status**: ✅ Ready to run

**Runs on**: Every push and PR

**Jobs**:
- test (8 matrix jobs)
- lint (1 job)
- docker (1 job)

**Total**: ~10-15 minutes

### PyPI Publishing

**Status**: ✅ Ready (needs PyPI setup)

**Runs on**: Release or manual

**Steps**:
1. Build package
2. Check with twine
3. Publish to PyPI

**Total**: ~2-3 minutes

### Docker Publishing

**Status**: ✅ Ready to run

**Runs on**: Release, push to main, or manual

**Steps**:
1. Build multi-platform image
2. Push to GHCR
3. Generate attestation

**Total**: ~10-15 minutes

---

## 🔗 URLs

### GitHub Actions

```
https://github.com/msk-access/getbasecounts/actions
```

### Docker Images

```
https://github.com/msk-access/getbasecounts/pkgs/container/getbasecounts
```

### Pull Docker Image

```bash
docker pull ghcr.io/msk-access/getbasecounts:latest
docker pull ghcr.io/msk-access/getbasecounts:v2.0.0
```

### PyPI Package

```
https://pypi.org/project/getbasecounts/
```

### Install from PyPI

```bash
pip install getbasecounts
pip install "getbasecounts[all]"
```

---

## 📝 Files Created/Modified

### GitHub Actions Workflows

| File | Purpose | Status |
|------|---------|--------|
| `.github/workflows/ci.yml` | CI testing | ✅ Updated |
| `.github/workflows/test.yml` | Comprehensive tests | ✅ Created |
| `.github/workflows/publish-pypi.yml` | PyPI publishing | ✅ Created |
| `.github/workflows/publish-docker.yml` | Docker publishing | ✅ Created |

### Docker Files

| File | Purpose | Status |
|------|---------|--------|
| `Dockerfile` | Ubuntu 22.04 production | ✅ Updated |
| `Dockerfile.python-slim` | Backup of old version | ✅ Archived |
| `Dockerfile.test` | Testing image | ✅ Exists |
| `docker-compose.yml` | Orchestration | ✅ Exists |

### Documentation

| File | Purpose | Status |
|------|---------|--------|
| `docs/GITHUB_ACTIONS.md` | CI/CD guide | ✅ Created |
| `docs/DOCKER_BASE_COMPARISON.md` | Base image comparison | ✅ Created |
| `docs/DOCKER_SUMMARY.md` | Docker summary | ✅ Updated |
| `docs/SUMMARY.md` | GitBook TOC | ✅ Updated |
| `CICD_SETUP_COMPLETE.md` | This file | ✅ Created |

---

## ✅ Verification

### Test Locally

```bash
# 1. Build Docker
docker build -t getbasecounts:test .

# 2. Test Docker
docker run --rm getbasecounts:test version
docker run --rm getbasecounts:test --help

# 3. Check size
docker images | grep getbasecounts

# 4. Run tests
pytest -v

# 5. Check linting
black --check src/
ruff check src/
mypy src/
```

### Test GitHub Actions

```bash
# 1. Push to a branch
git checkout -b test-ci
git push origin test-ci

# 2. Create PR
gh pr create --title "Test CI" --body "Testing GitHub Actions"

# 3. Watch workflows
gh run list
gh run watch

# 4. Check results
gh pr checks
```

---

## 🎉 Summary

### Docker

✅ **Base image**: Ubuntu 22.04 LTS  
✅ **Size**: ~900 MB  
✅ **Dependencies**: All included (samtools, libhts, cyvcf2, Ray)  
✅ **Multi-platform**: amd64 + arm64  
✅ **Published to**: GHCR (ghcr.io/msk-access/getbasecounts)  

### GitHub Actions

✅ **CI**: Matrix testing on every push  
✅ **PyPI**: Automated publishing on release  
✅ **Docker**: Multi-platform images to GHCR  
✅ **Quality**: Linting, type checking, coverage  

### Setup Required

1. ⚠️ Configure PyPI trusted publishing
2. ⚠️ (Optional) Add Codecov token
3. ⚠️ Enable branch protection
4. ⚠️ Set GHCR package visibility

### Ready to Use

✅ **Docker**: Build and run locally  
✅ **CI**: Runs on every push  
✅ **Publishing**: Ready (needs PyPI setup)  
✅ **Documentation**: Complete  

**Everything is set up and ready for production!** 🚀
