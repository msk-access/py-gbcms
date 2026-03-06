# py-gbcms (deprecated)

> ⚠️ **This package has been renamed to [`gbcms`](https://pypi.org/project/gbcms/).**

## Migration

```bash
pip uninstall py-gbcms
pip install gbcms
```

Update any `requirements.txt` or `pyproject.toml`:

```diff
-py-gbcms>=2.8.0
+gbcms>=3.0.0
```

This stub package (`py-gbcms==3.0.0`) exists solely to issue a deprecation
warning and pull in `gbcms>=3.0.0` as a dependency. It will not receive
further updates.
