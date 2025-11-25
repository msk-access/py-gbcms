# Documentation Development

## Local Preview

Preview documentation locally while editing:

```bash
# Install MkDocs (first time only)
pip install mkdocs-material

# Serve docs with live reload
mkdocs serve
```

Then open http://127.0.0.1:8000 in your browser. Changes to `docs/` files will auto-refresh.

## Documentation URLs

- **Production** (from `main`): https://msk-access.github.io/py-gbcms/
- **Staging** (from `develop`): https://msk-access.github.io/py-gbcms/gh-pages-staging/ (work in progress)
- **GitBook**: https://cmo-ci.gitbook.io/py-gbcms/

## Git-Flow Documentation Workflow

```bash
# 1. Work on feature branch
git checkout -b feature/update-docs

# 2. Edit docs and preview locally
mkdocs serve

# 3. Push to develop
git push origin feature/update-docs
# Create PR to develop

# 4. Merge to develop
# → Auto-deploys to staging URL

# 5. Release: merge develop → main  
# → Auto-deploys to production GitHub Pages
```

## Building Docs Manually

```bash
# Build static site
mkdocs build
# Output in site/ directory

# Deploy to GitHub Pages
mkdocs gh-deploy
```
