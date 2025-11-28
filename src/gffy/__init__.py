"""
gffy - GFF3 Genomic Statistics Calculator

A fast and efficient tool for computing comprehensive statistics from GFF3 genomic annotation files.
"""

__version__ = "0.0.1"
__author__ = "Emilio R"

from .compute_stats import compute_gff_stats
import os
import tempfile
import shutil
from urllib.request import urlopen
from contextlib import contextmanager
from collections import defaultdict, Counter
from .tools.constants import SKIP_FEATURES
from dataclasses import dataclass
import gzip

# -----------------------------
# Record class
# -----------------------------
@dataclass
class GFFRecord:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None
    strand: str
    phase: str
    attributes: dict

    @property
    def length(self) -> int:
        """Return feature length."""
        return self.end - self.start + 1


# -----------------------------
# Main Parser / Statistics Class
# -----------------------------
class GFFy:
    """
    A context-managed GFF parser and statistics calculator.
    Supports local files, URLs, and gzipped files.
    Computes counts, lengths (mean, median, min, max), spliced lengths,
    and allows custom filtering by feature, biotype, or other attributes.
    """

    # -----------------------
    # Initialization
    # -----------------------
    def __init__(self, source: str, skip_features: SKIP_FEATURES = SKIP_FEATURES, is_gzipped: bool = True):
        self.source = source
        self.skip_features = skip_features
        self.is_gzipped = is_gzipped
        self._temp_path = None  # For downloaded URLs
        self._file = None       # File handle
        self.records: list[GFFRecord] = []
        self.stats = None       # Placeholder for statistics object or cache

    # -----------------------
    # Context manager
    # -----------------------
    def __enter__(self):
        """Open file (local or URL) and parse GFF."""
        path = self._prepare_file(self.source)
        #self._file = self._open_stream(path, self.is_gzipped)
        #self._parse_stream()
        # Optionally, initialize cached statistics
        self.stats = {}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file and cleanup temporary files."""
        if self._file:
            self._file.close()
        if self._temp_path and os.path.exists(self._temp_path):
            os.remove(self._temp_path)
        return False  # propagate exceptions

    # -----------------------
    # Input handling / streaming
    # -----------------------
    @staticmethod
    def _is_url(s: str) -> bool:
        return isinstance(s, str) and s.startswith(("http://", "https://"))

    def _prepare_file(self, source: str) -> str:
        """Download URL to temp file if needed, or return local path."""
        if not self._is_url(source):
            return source
        fd, temp_path = tempfile.mkstemp(suffix=".gff" if not self.is_gzipped else ".gff.gz")
        os.close(fd)
        with urlopen(source) as r, open(temp_path, "wb") as f:
            shutil.copyfileobj(r, f)
        self._temp_path = temp_path
        return temp_path

    @staticmethod
    def _open_stream(path: str, is_gzipped: bool):
        """Open file or gzip file as text stream."""
        if is_gzipped:
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")

    # -----------------------
    # Parsing logic
    # -----------------------
    def _parse_stream(self):
        """Parse file line by line into GFFRecord objects."""
        for line in self._file:
            if line.startswith("#") or not line.strip():
                continue
            rec = self._parse_line(line)
            if rec:
                self.records.append(rec)

    def _parse_line(self, line: str) -> GFFRecord | None:
        """Parse a single GFF line into a GFFRecord."""
        # TODO: implement GFF parsing logic here
        pass

    # -----------------------
    # Core statistics functions
    # -----------------------
    def total_records(self) -> int:
        """Return total number of records."""
        return len(self.records)

    def feature_counts(self) -> dict:
        """Return counts per feature type."""
        return Counter(r.type for r in self.records)

    def biotype_counts(self) -> dict:
        """Return counts per biotype (from attributes)."""
        return Counter(r.attributes.get("biotype") for r in self.records)

    def mean_length(self, feature=None, biotype=None) -> float:
        """Return mean length of features, optionally filtered."""
        # TODO: implement filtering and mean computation
        pass

    def median_length(self, feature=None, biotype=None) -> float:
        """Return median feature length."""
        # TODO: implement median computation
        pass

    def min_length(self, feature=None, biotype=None) -> int:
        """Return minimum feature length."""
        # TODO: implement min computation
        pass

    def max_length(self, feature=None, biotype=None) -> int:
        """Return maximum feature length."""
        # TODO: implement max computation
        pass

    def spliced_length(self, transcript_id=None) -> dict:
        """Return spliced length per transcript (sum of exon lengths)."""
        # TODO: implement spliced length calculation
        pass

    def length_distribution(self, feature=None, biotype=None) -> list[int]:
        """Return list of lengths for plotting or analysis."""
        # TODO: implement
        pass

    # -----------------------
    # Custom filtering / query
    # -----------------------
    def filter_records(self, feature=None, biotype=None, seqid=None) -> list[GFFRecord]:
        """Return list of records matching filters."""
        # TODO: implement
        pass

    def count(self, feature=None, biotype=None, seqid=None) -> int:
        """Count records matching filters."""
        # TODO: implement
        pass

    def mean_filtered_length(self, feature=None, biotype=None, seqid=None) -> float:
        """Mean length of filtered records."""
        # TODO: implement
        pass

    def median_filtered_length(self, feature=None, biotype=None, seqid=None) -> float:
        """Median length of filtered records."""
        # TODO: implement
        pass

    def spliced_length_filtered(self, feature="transcript", biotype=None) -> dict:
        """Spliced length for filtered transcripts."""
        # TODO: implement
        pass

    def grouped_counts(self, attributes: list[str]) -> dict:
        """Counts grouped by multiple attributes."""
        # TODO: implement
        pass

    def grouped_mean_length(self, attributes: list[str]) -> dict:
        """Mean lengths grouped by multiple attributes."""
        # TODO: implement
        pass

    # -----------------------
    # QC / sanity checks
    # -----------------------
    def check_zero_length_features(self) -> list[GFFRecord]:
        """Return features with length = 0."""
        # TODO: implement
        pass

    def check_single_exon_genes(self) -> list[GFFRecord]:
        """Return genes with only one exon."""
        # TODO: implement
        pass

    def check_missing_cds(self) -> list[GFFRecord]:
        """Return transcripts missing CDS."""
        # TODO: implement
        pass

    def check_invalid_coordinates(self) -> list[GFFRecord]:
        """Return features with invalid start/end or strand."""
        # TODO: implement
        pass

    def summary_report(self) -> dict:
        """Return combined QC report."""
        # TODO: implement
        pass

    # -----------------------
    # Comparative / benchmarking
    # -----------------------
    def compare_feature_counts(self, other_parser) -> dict:
        """Compare counts per feature type with another parser."""
        # TODO: implement
        pass

    def compare_length_statistics(self, other_parser) -> dict:
        """Compare length statistics with another parser."""
        # TODO: implement
        pass

    def compare_spliced_lengths(self, other_parser) -> dict:
        """Compare spliced transcript lengths with another parser."""
        # TODO: implement
        pass

    def generate_comparison_table(self, other_parser, attributes=["feature", "biotype"]) -> dict:
        """Return combined table for benchmarking."""
        # TODO: implement
        pass

    # -----------------------
    # Convenience / utilities
    # -----------------------
    def to_dataframe(self):
        """Return records as pandas DataFrame."""
        # TODO: implement
        pass

    def histogram(self, feature=None, biotype=None, bins=50):
        """Return histogram data of lengths."""
        # TODO: implement
        pass

    def plot_length_distribution(self, feature=None, biotype=None):
        """Plot histogram of lengths."""
        # TODO: implement
        pass

    def export_summary(self, path, format="csv"):
        """Save summary statistics to file."""
        # TODO: implement
        pass

__all__ = ["GFFy", "compute_gff_stats"]
