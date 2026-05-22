"""
gffy - GFF3 Genomic Statistics Calculator

Fast, low-RAM statistics for GFF3 files. Output matches annotrieve GFFStats schema.

Usage
-----
**Command line** (after ``pip install gffy``):

    gffy annotation.gff3 --output stats.json
    gffy --from-file urls.txt -o stats.jsonl --workers 4
    gffy-convert annotation.gtf -o annotation.gff3
    python -m gffy annotation.gff3 --pretty

**Python library**:

    from gffy import compute_gff_stats, compute_bulk_stats, convert_gtf_to_gff3
    stats = compute_gff_stats("/path/to/annotation.gff3.gz")
    compute_bulk_stats(["a.gff3", "b.gff3"], "stats.jsonl", workers=2)
    convert_gtf_to_gff3("/path/to/annotation.gtf", "out.gff3")
    build_custom_annotation("/path/to/annotation.gff3")
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    __version__ = _pkg_version("gffy")
except PackageNotFoundError:
    __version__ = "0.0.0+local"

__author__ = "Emilio R"

from .bulk import (
    SourceListResult,
    compute_bulk_stats,
    read_source_list,
    warn_source_list_issues,
)
from .errors import ErrorInfo, classify_exception
from .cache import SourceInfo, build_sorted_gff, inspect_source, is_big
from .convert import convert_gtf_to_gff3, describe_convert_mode
from .annotrieve_export import build_custom_annotation, derive_custom_name
from .stats import compute_gff_stats, compute_gff_stats_from_lines, describe_compute_mode

__all__ = [
    "compute_gff_stats",
    "compute_gff_stats_from_lines",
    "compute_bulk_stats",
    "read_source_list",
    "SourceListResult",
    "warn_source_list_issues",
    "ErrorInfo",
    "classify_exception",
    "convert_gtf_to_gff3",
    "describe_compute_mode",
    "describe_convert_mode",
    "inspect_source",
    "is_big",
    "build_sorted_gff",
    "SourceInfo",
    "build_custom_annotation",
    "derive_custom_name",
    "__version__",
]
