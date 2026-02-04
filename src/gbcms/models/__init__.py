"""
Data models for gbcms.

Provides Pydantic models for variants, configuration, and core data structures.
"""

from .core import GbcmsConfig, GenomicInterval, OutputFormat, Variant, VariantType

__all__ = [
    "GbcmsConfig",
    "GenomicInterval",
    "OutputFormat",
    "Variant",
    "VariantType",
]
