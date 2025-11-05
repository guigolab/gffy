"""
gffy - GFF3 Genomic Statistics Calculator

A fast and efficient tool for computing comprehensive statistics from GFF3 genomic annotation files.
"""

__version__ = "0.1.0"
__author__ = "Emilio R"

from .compute_stats import compute_gff_stats

__all__ = ["compute_gff_stats"]

