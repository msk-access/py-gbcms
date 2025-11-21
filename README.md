# py-gbcms Documentation

This branch contains the documentation for py-gbcms, built with GitBook.

## 📖 Viewing the Documentation

**Live documentation:** https://msk-access.github.io/py-gbcms/

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

3. **Test locally:**
   ```bash
   cd docs
   npm install -g gitbook-cli
   gitbook install
   gitbook serve
   # Visit http://localhost:4000
   ```

4. **Commit and push:**
   ```bash
   git add docs/
   git commit -m "docs: Update quick start guide"
   git push origin gh-pages
   ```

5. **Automatic deployment:**
   - GitHub Actions will automatically build and deploy
   - Changes appear at https://msk-access.github.io/py-gbcms/ in ~2-3 minutes

### Adding New Pages

1. Create new markdown file in `docs/`:
   ```bash
   echo "# My New Page" > docs/my-new-page.md
   ```

2. Add to `docs/SUMMARY.md`:
   ```markdown
   * [My New Page](my-new-page.md)
   ```

3. Commit and push

## 📁 Structure

```
gh-pages/
├── .github/
│   └── workflows/
│       └── deploy-docs.yml    # Auto-deployment workflow
├── docs/                      # Documentation source (EDIT HERE)
│   ├── README.md
│   ├── SUMMARY.md             # Table of contents
│   ├── book.json              # GitBook config
│   ├── quick-start.md
│   ├── CLI_FEATURES.md
│   └── ...
├── .gitbook.yaml              # GitBook settings
└── [Built HTML files]         # Auto-generated (DO NOT EDIT)
```

## 🔧 GitBook Configuration

- **Main config:** `docs/book.json`
- **Structure:** `docs/SUMMARY.md`
- **Plugins:** Defined in `book.json`

## 📝 Documentation Guidelines

- Use clear, concise language
- Include code examples
- Add screenshots where helpful
- Keep navigation (SUMMARY.md) organized
- Test locally before pushing

## 🚀 Deployment

Deployment is automatic via GitHub Actions:
- **Trigger:** Push to `gh-pages` branch affecting `docs/`
- **Process:** GitBook build → Deploy to gh-pages root
- **Time:** ~2-3 minutes
- **URL:** https://msk-access.github.io/py-gbcms/

## 📮 Questions?

For documentation issues, open an issue on the main repository:
https://github.com/msk-access/py-gbcms/issues
