# GitHub Actions CI/CD

Complete guide to the GitHub Actions workflows for GetBaseCounts.

## Overview

GetBaseCounts uses GitHub Actions for:
- ✅ Continuous Integration (CI) - Testing on every push
- ✅ Publishing to PyPI - Automated package releases
- ✅ Publishing to GHCR - Docker image distribution
- ✅ Multi-platform testing - Linux and macOS
- ✅ Code quality checks - Linting and type checking

---

## Workflows

### 1. CI Workflow (`ci.yml`)

**Trigger**: Push or PR to `main` or `develop` branches

**Jobs**:
1. **test** - Run tests on multiple OS and Python versions
2. **lint** - Code quality checks (black, ruff, mypy)
3. **docker** - Build and test Docker images

**Matrix Testing**:
- OS: Ubuntu, macOS
- Python: 3.9, 3.10, 3.11, 3.12

**What it does**:
```yaml
✅ Install system dependencies (samtools, libhts)
✅ Install Python dependencies with all extras
✅ Run pytest with coverage
✅ Upload coverage to Codecov
✅ Run linters (black, ruff, mypy)
✅ Build Docker images
✅ Test Docker images
```

**Badge**:
```markdown
![CI](https://github.com/msk-access/getbasecounts/workflows/CI/badge.svg)
```

---

### 2. Test Workflow (`test.yml`)

**Trigger**: Push, PR, or manual dispatch

**Jobs**:
1. **test** - Comprehensive testing
2. **lint** - Code quality
3. **docker-test** - Docker validation

**Features**:
- Matrix testing across OS and Python versions
- Coverage reporting to Codecov
- Docker build caching
- Multi-platform Docker builds

---

### 3. Publish to PyPI (`publish-pypi.yml`)

**Trigger**: 
- Push tags matching `[0-9]+.[0-9]+.[0-9]+` (e.g., 1.0.0, 2.3.4)
- Release published
- Manual workflow dispatch

**What it does**:
```yaml
✅ Build Python package (wheel + sdist)
✅ Verify package with twine
✅ Publish to PyPI (only on tag push)
✅ Uses trusted publishing (no tokens needed)
✅ Skip existing versions
✅ Verbose output with package hash
```

**Automatic on Tag Push**:
```bash
# Create and push a version tag (without 'v' prefix)
git tag 2.0.0
git push origin 2.0.0
# Workflow automatically publishes to PyPI
```

**Or via Release**:
```bash
# Create a release on GitHub
gh release create 2.0.0 --title "Version 2.0.0" --notes "Release notes"
# Workflow automatically publishes to PyPI
```

**Manual Trigger**:
```bash
# Go to Actions tab → Publish Python Package → Run workflow
# Optional: Specify version in input
```

**Setup Required**:
1. Configure PyPI trusted publishing:
   - Go to https://pypi.org/manage/account/publishing/
   - Add GitHub repository: `msk-access/getbasecounts`
   - Workflow: `publish-pypi.yml`
   - Environment: (leave blank)

---

### 4. Build, Test, and Push Docker Image (`publish-docker.yml`)

**Trigger**:
- Push any tag
- Pull requests to `main`, `master`, or `develop`
- Manual workflow dispatch

**Jobs**:
1. **build** - Build Docker image (runs on all triggers)
2. **push** - Push to GHCR (only on tag push)

**What it does**:
```yaml
✅ Build Docker image (Ubuntu 22.04 base)
✅ Test build on PRs (no push)
✅ Push to GitHub Container Registry (only on tags)
✅ Automatic tagging (latest, tag name)
✅ Uses GITHUB_TOKEN (no secrets needed)
```

**Image Tags**:
- `latest` - Latest tagged version
- `2.0.0` - Specific tag name

**Pull Image**:
```bash
docker pull ghcr.io/msk-access/getbasecounts:latest
docker pull ghcr.io/msk-access/getbasecounts:2.0.0
```

**How to Publish**:
```bash
# Create and push a tag
git tag 2.0.0
git push origin 2.0.0

# Workflow automatically:
# 1. Builds Docker image
# 2. Pushes to ghcr.io/msk-access/getbasecounts:2.0.0
# 3. Tags as latest
```

**Setup Required**:
- None! Uses automatic `GITHUB_TOKEN`

---

## Setup Instructions

### 1. Enable GitHub Actions

1. Go to repository Settings → Actions → General
2. Set "Actions permissions" to "Allow all actions"
3. Set "Workflow permissions" to "Read and write permissions"
4. Save

### 2. Configure PyPI Publishing

**Option A: Trusted Publishing (Recommended)**

1. Go to https://pypi.org/manage/account/publishing/
2. Click "Add a new publisher"
3. Fill in:
   - PyPI Project Name: `getbasecounts`
   - Owner: `msk-access`
   - Repository: `getbasecounts`
   - Workflow: `publish-pypi.yml`
   - Environment: (leave blank)
4. Save

**Option B: API Token**

1. Generate PyPI API token
2. Add to GitHub Secrets:
   - Name: `PYPI_API_TOKEN`
   - Value: `pypi-...`
3. Update workflow to use token

### 3. Configure Codecov (Optional)

1. Go to https://codecov.io/
2. Add repository
3. Get upload token
4. Add to GitHub Secrets:
   - Name: `CODECOV_TOKEN`
   - Value: `...`

### 4. Configure GitHub Container Registry

**Automatic** - No setup needed!

GitHub Actions automatically has permission to push to GHCR.

---

## Usage

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
```

### Manual Workflow Triggers

**Publish to Test PyPI**:
```bash
# Via GitHub UI
1. Go to Actions tab
2. Select "Publish to PyPI"
3. Click "Run workflow"
4. Check "Publish to Test PyPI"
5. Click "Run workflow"
```

**Publish Docker Image**:
```bash
# Via GitHub UI
1. Go to Actions tab
2. Select "Publish Docker Image"
3. Click "Run workflow"
4. Select branch
5. Click "Run workflow"
```

### Creating a Release

```bash
# Using GitHub CLI
gh release create v2.0.0 \
  --title "Version 2.0.0" \
  --notes "Release notes here"

# This triggers:
# 1. publish-pypi.yml → Publishes to PyPI
# 2. publish-docker.yml → Publishes Docker image
```

---

## Workflow Details

### CI Workflow Matrix

```yaml
OS: [ubuntu-latest, macos-latest]
Python: [3.9, 3.10, 3.11, 3.12]
Total: 8 jobs
```

**Why this matrix?**
- Ubuntu: Most common deployment target
- macOS: Developer machines
- Python 3.9+: Current supported versions

### System Dependencies

**Ubuntu**:
```bash
build-essential
zlib1g-dev
libbz2-dev
liblzma-dev
libcurl4-openssl-dev
libssl-dev
libhts-dev  # For cyvcf2
samtools
```

**macOS**:
```bash
htslib  # For cyvcf2
samtools
```

### Python Dependencies

```bash
# Core + all optional features
uv pip install -e ".[dev,all]"

# Includes:
# - Core dependencies (pysam, numpy, etc.)
# - cyvcf2 (fast VCF parsing)
# - Ray (distributed computing)
# - Dev tools (pytest, black, ruff, mypy)
```

---

## Docker Build Process

### Multi-Stage Build

```dockerfile
Stage 1 (builder):
- Ubuntu 22.04
- Install build dependencies
- Install Python 3.11
- Build package with all features
- Size: ~1.5 GB (discarded)

Stage 2 (final):
- Ubuntu 22.04
- Install runtime dependencies only
- Copy built packages from stage 1
- Size: ~900 MB
```

### Multi-Platform Support

```yaml
platforms: linux/amd64,linux/arm64
```

**Why?**
- amd64: Intel/AMD processors (most servers)
- arm64: Apple Silicon, AWS Graviton

### Caching Strategy

```yaml
cache-from: type=gha
cache-to: type=gha,mode=max
```

**Benefits**:
- Faster builds (reuse layers)
- Lower bandwidth usage
- Faster CI/CD pipelines

---

## Troubleshooting

### Issue: Tests Fail on macOS

**Problem**: Missing system dependencies

**Solution**:
```bash
brew install htslib samtools
```

### Issue: Docker Build Fails

**Problem**: Out of disk space

**Solution**: Clean up Docker
```bash
docker system prune -a
```

### Issue: PyPI Publishing Fails

**Problem**: Package already exists

**Solution**: Increment version in `pyproject.toml`

### Issue: GHCR Push Fails

**Problem**: Permission denied

**Solution**: Check workflow permissions
1. Settings → Actions → General
2. Set "Workflow permissions" to "Read and write"

### Issue: Coverage Upload Fails

**Problem**: Codecov token missing

**Solution**: Add `CODECOV_TOKEN` to secrets or use tokenless upload

---

## Best Practices

### 1. Version Tagging

```bash
# Use semantic versioning
git tag v2.0.0
git push origin v2.0.0

# This triggers Docker build with proper tags
```

### 2. Pre-Release Testing

```bash
# Test on Test PyPI first
1. Manual trigger workflow
2. Check "Publish to Test PyPI"
3. Test installation: pip install -i https://test.pypi.org/simple/ getbasecounts
4. If OK, create release for real PyPI
```

### 3. Docker Image Verification

```bash
# After push, verify image
docker pull ghcr.io/msk-access/getbasecounts:latest
docker run --rm ghcr.io/msk-access/getbasecounts:latest version
```

### 4. Branch Protection

Enable branch protection for `main`:
- Require status checks to pass
- Require CI workflow to pass
- Require code review

---

## Monitoring

### Check Workflow Status

```bash
# Via GitHub CLI
gh run list
gh run view <run-id>
gh run watch <run-id>
```

### View Logs

```bash
# Via GitHub CLI
gh run view <run-id> --log

# Or via web UI
https://github.com/msk-access/getbasecounts/actions
```

### Coverage Reports

View on Codecov:
```
https://codecov.io/gh/msk-access/getbasecounts
```

### Docker Images

View on GHCR:
```
https://github.com/msk-access/getbasecounts/pkgs/container/getbasecounts
```

---

## Workflow Files

| File | Purpose | Trigger |
|------|---------|---------|
| `ci.yml` | Continuous Integration | Push, PR |
| `test.yml` | Comprehensive Testing | Push, PR, Manual |
| `publish-pypi.yml` | PyPI Publishing | Release, Manual |
| `publish-docker.yml` | Docker Publishing | Release, Push to main, Tags |

---

## Security

### Secrets Used

| Secret | Purpose | Required |
|--------|---------|----------|
| `GITHUB_TOKEN` | Automatic | ✅ Auto-provided |
| `CODECOV_TOKEN` | Coverage upload | ⚠️ Optional |
| `PYPI_API_TOKEN` | PyPI (if not using trusted publishing) | ⚠️ Optional |

### Permissions

```yaml
contents: read      # Read repository
packages: write     # Push to GHCR
id-token: write     # Trusted publishing
```

### Build Attestation

Docker images include build provenance:
- Source repository
- Commit SHA
- Workflow run
- Build timestamp

---

## Summary

### Workflows

✅ **CI** - Test on every push  
✅ **Test** - Comprehensive testing  
✅ **PyPI** - Automated package publishing  
✅ **Docker** - Multi-platform image publishing  

### Features

✅ **Matrix testing** - Multiple OS and Python versions  
✅ **Code quality** - Linting and type checking  
✅ **Coverage** - Codecov integration  
✅ **Caching** - Fast builds  
✅ **Security** - Trusted publishing, attestation  
✅ **Multi-platform** - amd64 and arm64  

### Setup

1. Enable GitHub Actions
2. Configure PyPI trusted publishing
3. (Optional) Add Codecov token
4. Create releases to trigger publishing

**All workflows are ready to use!** 🚀
