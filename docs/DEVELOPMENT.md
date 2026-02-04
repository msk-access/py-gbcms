# Developer Guide

Guide for contributing to py-gbcms.

---

## Setup

```bash
# Clone
git clone https://github.com/msk-access/py-gbcms.git
cd py-gbcms

# Virtual environment
python -m venv .venv
source .venv/bin/activate

# Install (builds Rust extension)
maturin develop --release

# Verify
gbcms --version
```

---

## Project Structure

```mermaid
flowchart LR
    subgraph Python["src/gbcms/"]
        CLI[cli.py] --> Pipeline[pipeline.py]
        Pipeline --> IO[io/]
        Pipeline --> Models[models/]
    end
    
    subgraph Rust["rust/src/"]
        Lib[lib.rs] --> Count[counting.rs]
        Count --> Stats[stats.rs]
    end
    
    Pipeline --> Rust
```

---

## Build Commands

```bash
# Development (fast)
maturin develop

# Release (optimized)
maturin develop --release

# Build wheel
maturin build --release --out dist
```

---

## Code Standards

### Python

| Standard | Requirement |
|:---------|:------------|
| Type hints | All public functions |
| Docstrings | Google style |
| Exports | `__all__` in every module |
| Logging | Use `logging`, not `print()` |
| Config | Pydantic models |

### Rust

| Standard | Requirement |
|:---------|:------------|
| Docs | `///` on public items |
| Errors | `anyhow::Result` |
| Logging | `log` crate |

---

## Git Workflow (git-flow)

```mermaid
gitGraph
    commit id: "main"
    branch develop
    commit id: "develop"
    branch feature/new-thing
    commit id: "work"
    checkout develop
    merge feature/new-thing
    branch release/2.3.0
    commit id: "bump"
    checkout main
    merge release/2.3.0 tag: "2.3.0"
    checkout develop
    merge release/2.3.0
```

| Branch | Purpose |
|:-------|:--------|
| `main` | Production releases |
| `develop` | Integration |
| `feature/*` | New features |
| `release/*` | Release candidates |
| `hotfix/*` | Production fixes |

---

## Quality Checklist

Before committing:

- [ ] `make lint` passes
- [ ] `pytest` passes
- [ ] Type hints complete
- [ ] Docstrings added
- [ ] No dead code

---

## Environment Variables

| Variable | Default | Description |
|:---------|:--------|:------------|
| `GBCMS_LOG_LEVEL` | INFO | Logging level |
| `RUST_LOG` | — | Rust logging |

```bash
GBCMS_LOG_LEVEL=DEBUG RUST_LOG=debug gbcms run ...
```
