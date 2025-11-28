"""
gffy - GFF3 Genomic Statistics Calculator

A fast and efficient tool for computing comprehensive statistics from GFF3 genomic annotation files.
"""

__version__ = "0.0.1"
__author__ = "Emilio R"

from .compute_stats import compute_gff_stats
from .tools.classes import Transcript, Gene
from .tools.constants import SKIP_FEATURES, FEATURE_CODES
from .tools import helpers
from collections import Counter
import os
import sys


class GFFy:
    """
    A context-managed GFF parser and statistics calculator.
    Supports local files, URLs, and gzipped files.
    """

    def __init__(
        self, 
        source: str, 
        skip_features: frozenset = SKIP_FEATURES,
        is_gzipped: bool | str = 'auto'
    ):
        self.source = source
        self.skip_features = skip_features
        self.is_gzipped = is_gzipped if is_gzipped != 'auto' else helpers.is_gzipped(source)
        self._temp_path = None
        self._file = None
        self.transcripts: dict[str, Transcript] = {}
        self.biotypes: Counter[str] = Counter()
        self.feature_types: Counter[str] = Counter()
        self.genes: dict[str, Gene] = {}
        self.stats = None

    def __enter__(self):
        """Open file (local or URL) and parse GFF."""
        path = self._prepare_file()
        self._file = helpers.open_gff_stream(path, self.is_gzipped)
        self._parse_transcripts()
        self.stats = {}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close file and cleanup temporary files."""
        if self._file:
            self._file.close()
        if self._temp_path and os.path.exists(self._temp_path):
            os.remove(self._temp_path)
        return False

    def _prepare_file(self) -> str:
        """Download URL to temp file if needed, or return local path."""
        if not helpers.is_url(self.source):
            return self.source
        self._temp_path = helpers.download_url_to_temp(self.source, self.is_gzipped)
        return self._temp_path

    def _parse_transcripts(self):
        """Parse GFF file and populate transcripts dictionary."""
        # Cache frequently accessed methods
        feat_code_get = FEATURE_CODES.get
        transcripts = self.transcripts
        intern = sys.intern
        
        for line in helpers.stream_gff_lines(self._file):
            # Split line into columns
            cols = line.split("\t", 8)
            if len(cols) < 9:
                continue
            
            # Quick check: must have Parent= in attributes column
            attrs = cols[8]
            if "Parent=" not in attrs:
                continue
            
            # Parse basic feature info
            feature_type = intern(cols[2])
            start = int(cols[3])
            end = int(cols[4])
            length = end - start + 1
            
            # Get feature code and parse attributes
            feature_code = feat_code_get(feature_type, 0)
            feature_id, parent_id, biotype = helpers.parse_gff_attributes_fast(attrs, feature_code)
            
            if not parent_id:
                continue
            
            # Process transcript feature (mRNA, lncRNA, etc.)
            if feature_code == 0 and feature_id:
                helpers.process_transcript_feature(
                    feature_id, parent_id, feature_type, biotype, length, transcripts
                )
            
            # Process exon or CDS
            elif feature_code != 0:
                helpers.process_exon_or_cds(parent_id, feature_code, length, transcripts)

__all__ = ["GFFy", "compute_gff_stats"]
