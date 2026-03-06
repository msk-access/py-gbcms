"""
py-gbcms is deprecated.

The package has been renamed to `gbcms`. Please update your dependencies:

    pip install gbcms

    # requirements.txt / pyproject.toml:
    # gbcms>=3.0.0  (replace py-gbcms)

This shim package re-exports everything from gbcms so existing code
continues to work, but the deprecation warning will appear on import.
"""

import warnings

warnings.warn(
    "py-gbcms has been renamed to gbcms. "
    "Please update your dependency: `pip install gbcms`. "
    "This shim will be removed in a future release.",
    DeprecationWarning,
    stacklevel=2,
)

# Re-export everything from the real package so existing imports keep working
from gbcms import *  # noqa: F401, F403
from gbcms import __version__  # noqa: F401
