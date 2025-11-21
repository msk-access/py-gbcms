# py-gbcms Documentation

This branch contains the documentation for py-gbcms.

## �� Viewing the Documentation

Documentation is available at **two locations**:

1. **GitBook.com** (Primary): https://cmo-ci.gitbook.io/py-gbcms/
2. **GitHub Pages** (Mirror): https://msk-access.github.io/py-gbcms/

Both sites are built from the same markdown files in the `docs/` folder.

## ✏️ Contributing to Documentation

All documentation lives in this `gh-pages` branch in the `docs/` folder.

### Making Changes

1. **Clone and checkout:**
   ```bash
   git clone https://github.com/msk-access/py-gbcms.git
   cd py-gbcms
   git checkout gh-pages
   ```

2. **Edit markdown files:**
   ```bash
   vim docs/quick-start.md
   vim docs/CLI_FEATURES.md
   # etc.
   ```

3. **Preview locally (GitBook - Optional):**
   ```bash
   cd docs
   npm install -g gitbook-cli
   gitbook serve
   # Visit http://localhost:4000
   ```

4. **Preview locally (MkDocs - Recommended):**
   ```bash
   pip install mkdocs-material
   mkdocs serve
   # Visit http://127.0.0.1:8000
   ```

5. **Commit and push:**
   ```bash
   git add docs/
   git commit -m "docs: Update quick start guide"
   git push origin gh-pages
   ```

6. **Automatic deployment:**
   - **GitBook.com**: Updates automatically from gh-pages branch
   - **GitHub Pages**: GitHub Actions builds with MkDocs (~2-3 minutes)

### Adding New Pages

1. Create new markdown file in `docs/`:
   ```bash
   echo "# My New Page" > docs/my-new-page.md
   ```

2. Add to `docs/SUMMARY.md` (for GitBook):
   ```markdown
   * [My New Page](my-new-page.md)
   ```

3. Add to `mkdocs.yml` (for GitHub Pages):
   ```yaml
   nav:
     - My New Page: my-new-page.md
   ```

4. Commit and push

## 📁 Structure

```
gh-pages/
├── .github/
│   └── workflows/
│       └── deploy-docs.yml    # MkDocs → GitHub Pages
├── docs/                      # Documentation source (EDIT HERE)
│   ├── README.md
│   ├── SUMMARY.md             # GitBook table of contents
│   ├── book.json              # GitBook config
│   ├── CONTRIBUTING.md
│   ├── CHANGELOG.md
│   ├── quick-start.md
│   ├── CLI_FEATURES.md
│   └── ...
├── .gitbook.yaml              # GitBook.com settings
├── mkdocs.yml                 # MkDocs settings (GitHub Pages)
└── README.md                  # This file
```

## 🔧 Configuration Files

### For GitBook.com
- **Main config:** `docs/book.json`
- **Structure:** `docs/SUMMARY.md`
- **Settings:** `.gitbook.yaml`

### For GitHub Pages (MkDocs)
- **Main config:** `mkdocs.yml`
- **Theme:** Material for MkDocs
- **Deployment:** `.github/workflows/deploy-docs.yml`

## 📝 Documentation Guidelines

- Use clear, concise language
- Include code examples
- Add screenshots where helpful
- Keep navigation (SUMMARY.md & mkdocs.yml) in sync
- Test locally before pushing

## 🚀 Deployment

### GitBook.com
- **Trigger:** Automatic sync from gh-pages branch
- **URL:** https://cmo-ci.gitbook.io/py-gbcms/

### GitHub Pages
- **Trigger:** Push to `gh-pages` branch affecting `docs/` or `mkdocs.yml`
- **Process:** MkDocs build → Deploy to gh-pages (in `gh-deploy` commit)
- **Time:** ~2-3 minutes
- **URL:** https://msk-access.github.io/py-gbcms/

## 📮 Questions?

For documentation issues, open an issue on the main repository:
https://github.com/msk-access/py-gbcms/issues
