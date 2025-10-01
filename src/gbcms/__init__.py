"""
GetBaseCounts - Calculate base counts in multiple BAM files for variants

A high-performance Python reimplementation of GetBaseCountsMultiSample.
"""

__version__ = "2.0.0"
__author__ = "MSK-ACCESS Team"

from .config import Config, CountType
from .variant import VariantEntry

__all__ = ["Config", "CountType", "VariantEntry", "__version__"]
