# Contributing to py-gbcms

Thank you for your interest in contributing to py-gbcms! This guide covers both code contributions and documentation updates.

## Table of Contents

- [Code Contributions](#code-contributions)
- [Documentation Contributions](#documentation-contributions)
- [Reporting Issues](#reporting-issues)

---

## Code Contributions

### Development Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/msk-access/py-gbcms.git
   cd py-gbcms
   ```

2. **Install uv** (if not already installed)
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

3. **Install development dependencies**
   ```bash
   uv pip install -e ".[dev]"
   ```

### Development Workflow

#### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=gbcms --cov-report=html

# Run specific test file
pytest tests/test_counter.py

# Run with verbose output
pytest -v
```

#### Code Quality

```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/

# Fix linting issues automatically
ruff check --fix src/ tests/

# Type checking
mypy src/
```

#### Building and Testing Docker

```bash
# Build Docker image
docker build -t gbcms:latest .
```

### Code Style

- Follow PEP 8 style guidelines
- Use type hints for all function signatures
- Write docstrings for all public functions and classes
- Keep functions focused and single-purpose
- Maximum line length: 100 characters

### Testing Guidelines

- Write tests for all new features
- Maintain or improve code coverage
- Use descriptive test names
- Use fixtures for common test setup
- Test edge cases and error conditions

### Pull Request Process

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Write code following the style guidelines
   - Add tests for new functionality
   - Update documentation as needed

3. **Run tests and linters**
   ```bash
   make test
   make lint
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Description of your changes"
   ```

5. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request**
   - Provide a clear description of the changes
   - Reference any related issues
   - Ensure all CI checks pass

---

## Documentation Contributions

### Local Preview

Preview documentation locally while editing:

```bash
# Install MkDocs (first time only)
pip install mkdocs-material

# Serve docs with live reload
mkdocs serve
```

Then open http://127.0.0.1:8000 in your browser. Changes to `docs/` files will auto-refresh.

### Documentation URLs

- **Production** (from `main`): https://msk-access.github.io/py-gbcms/
- **Staging** (from `develop`): https://msk-access.github.io/py-gbcms/gh-pages-staging/
- **GitBook**: https://cmo-ci.gitbook.io/py-gbcms/

### Git-Flow Documentation Workflow

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

### Building Docs Manually

```bash
# Build static site
mkdocs build
# Output in site/ directory

# Deploy to GitHub Pages (maintainers only)
mkdocs gh-deploy
```

---

## Reporting Issues

When reporting issues, please include:

- Python version
- Operating system
- Steps to reproduce
- Expected behavior
- Actual behavior
- Error messages or logs

## Feature Requests

We welcome feature requests! Please:

- Check if the feature already exists
- Provide a clear use case
- Describe the expected behavior
- Consider submitting a PR if you can implement it

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Help others learn and grow

## Questions?

Feel free to open an issue for questions or reach out to the maintainers.

Thank you for contributing!
