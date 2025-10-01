# GitHub Workflows Updated ✅

Summary of all GitHub Actions workflow updates for GetBaseCounts CI/CD pipelines.

---

## ✅ What Was Updated

### 1. PyPI Publishing Workflow

**File**: `.github/workflows/publish-pypi.yml`

**Changes**:
- ✅ Updated to follow best practices for Python package publishing
- ✅ Triggers on tag push matching `[0-9]+.[0-9]+.[0-9]+` (e.g., `2.0.0`)
- ✅ Also triggers on release or manual dispatch
- ✅ Only publishes when tags are pushed
- ✅ Uses pip cache for faster builds
- ✅ Fetches full git history
- ✅ More verbose output with package hash
- ✅ Shows package info after build

**Key Features**:
- Tag format: `2.0.0` (NO 'v' prefix)
- Publishing only happens on tag push
- Trusted publishing (no API tokens needed)

### 2. Docker Publishing Workflow

**File**: `.github/workflows/publish-docker.yml`

**Changes**:
- ✅ Updated for GetBaseCounts Docker publishing
- ✅ Renamed to "Build, Test, and Push Docker Image"
- ✅ Split into two jobs: `build` and `push`
- ✅ Build job runs on all triggers (PRs, tags, manual)
- ✅ Push job only runs on tag push
- ✅ Simplified tagging (latest + tag name only)
- ✅ Uses GITHUB_TOKEN (no custom secrets needed)

**Key Features**:
- Tests build on every PR
- Only pushes on tag push
- Ubuntu 22.04 base image
- All dependencies included (samtools, libhts, cyvcf2, Ray)

---

## 🚀 How to Use

### Publishing a New Version

**Step 1: Update Version**
```bash
# Edit pyproject.toml
# Change version = "2.0.0"
git add pyproject.toml
git commit -m "Bump version to 2.0.0"
git push
```

**Step 2: Create and Push Tag**
```bash
# Create tag WITHOUT 'v' prefix
git tag 2.0.0
git push origin 2.0.0
```

**Step 3: Automatic Publishing**
```
This automatically triggers:
1. publish-pypi.yml → Publishes to PyPI
2. publish-docker.yml → Builds and pushes Docker image to GHCR

Results:
- PyPI: pip install getbasecounts==2.0.0
- GHCR: docker pull ghcr.io/msk-access/getbasecounts:2.0.0
- GHCR: docker pull ghcr.io/msk-access/getbasecounts:latest
```

---

## 📊 Workflow Triggers

### publish-pypi.yml

| Trigger | Action |
|---------|--------|
| Push tag `2.0.0` | ✅ Build + Publish to PyPI |
| Release published | ✅ Build + Publish to PyPI |
| Manual dispatch | ✅ Build only (no publish) |

### publish-docker.yml

|---------|-----------|----------|
| Push any tag | ✅ Runs | ✅ Runs (pushes to GHCR) |
| Pull request | ✅ Runs | ❌ Skipped |
| Manual dispatch | ✅ Runs | ❌ Skipped |

## 🎯 Workflow Files

| File | Purpose | Status |
|------|---------|--------|
| `.github/workflows/ci.yml` | CI testing | ✅ Updated |
| `.github/workflows/test.yml` | Comprehensive tests | ✅ Updated |
| `.github/workflows/publish-pypi.yml` | PyPI publishing | ✅ Updated |
| `.github/workflows/publish-docker.yml` | Docker publishing | ✅ Updated |

---

## ✅ Verification Checklist

- [ ] PyPI trusted publishing configured
- [ ] Repository has write permissions for Actions
- [ ] Version in `pyproject.toml` is correct
- [ ] All tests pass locally
- [ ] Docker builds locally

**All workflows are ready for production!** 🚀
