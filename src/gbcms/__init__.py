"""
gbcms (Get Base Counts Multi-Sample) - A tool for counting bases at variant positions.

This package provides a command-line interface and Python API for genotyping
variants in BAM files using a high-performance Rust counting engine.

Example usage:
    $ gbcms run -v variants.vcf -b sample.bam -f reference.fa -o output/
"""

__version__ = "2.7.0"

from .models.core import GbcmsConfig, OutputFormat, Variant, VariantType
from .pipeline import Pipeline

__all__ = [
    "__version__",
    "GbcmsConfig",
    "OutputFormat",
    "Pipeline",
    "Variant",
    "VariantType",
]
