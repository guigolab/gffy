from collections import defaultdict
from typing import Optional, IO
from .classes import Gene, Transcript, OrphanFeature
from .constants import COMMON_TRANSCRIPT_TYPES, COMMON_GENE_TYPES, FEATURE_CODES
import sys
import requests
import gzip
import os
import tempfile
import shutil
from urllib.request import urlopen
from contextlib import contextmanager
from array import array

def init_orphan_feature(
    feature_id: str,
    feature_type: str,
    length: int,
    parent_id: Optional[str],
    biotype: Optional[str],
) -> OrphanFeature:
    return OrphanFeature(feature_id, feature_type, length, parent_id, biotype)

def categorize_roots(roots: dict) -> dict:
    root_to_category = {}
    for root_id, info in roots.items():
        feature_type = info.feature_type
        biotype = info.biotype or ""
        if feature_type == "pseudogene":
            root_to_category[root_id] = "pseudogene"
        elif info.has_cds or "protein_coding" in biotype.lower():
            root_to_category[root_id] = "coding"
        elif info.has_exon:
            # Classify non-coding genes based on length and exon count
            gene_length = info.length
            
            # Long non-coding: >200bp AND multiple exons
            # Short non-coding: <=200bp OR single exon
            if gene_length > 200 or info.has_multiple_exons:
                root_to_category[root_id] = "long_non_coding"
            else:
                root_to_category[root_id] = "short_non_coding"
        else:
            root_to_category[root_id] = None
    return root_to_category

    
def insert_exon_sorted(transcript: Transcript, start: int, end: int) -> None:
    """Insert an exon interval while keeping the flattened array sorted by start."""
    insert_idx = transcript.exon_count
    coords = transcript.exons_lengths

    if insert_idx == 0 or start >= coords[(insert_idx - 1) * 2]:
        coords.append(start)
        coords.append(end)
    else:
        lo = 0
        hi = insert_idx
        while lo < hi:
            mid = (lo + hi) // 2
            mid_start = coords[mid * 2]
            if start >= mid_start:
                lo = mid + 1
            else:
                hi = mid
        insert_pos = lo * 2
        coords.insert(insert_pos, start)
        coords.insert(insert_pos + 1, end)

    transcript.exon_count += 1


def init_category_totals() -> dict:
    return {
        "coding": {
            "exons": 0,
            "introns": 0,
            "cds": 0,
            "exon_len_sum": 0,
            "intron_len_sum": 0,
            "cds_len_sum": 0,
            "cds_transcripts": 0,
        },
        "long_non_coding": {
            "exons": 0,
            "introns": 0,
            "cds": 0,
            "exon_len_sum": 0,
            "intron_len_sum": 0,
            "cds_len_sum": 0,
            "cds_transcripts": 0,
        },
        "short_non_coding": {
            "exons": 0,
            "introns": 0,
            "cds": 0,
            "exon_len_sum": 0,
            "intron_len_sum": 0,
            "cds_len_sum": 0,
            "cds_transcripts": 0,
        },
        "pseudogene": {
            "exons": 0,
            "introns": 0,
            "cds": 0,
            "exon_len_sum": 0,
            "intron_len_sum": 0,
            "cds_len_sum": 0,
            "cds_transcripts": 0,
        },
    }


def parse_gff_line_fast(cols: list) -> tuple:
    feature_type = sys.intern(cols[2])
    start = int(cols[3])
    end = int(cols[4])
    length = end - start + 1
    attr_str = cols[8].strip("\n")
    parent_id = None
    biotype = None
    feature_id = None
    for attr in attr_str.split(";"):
        #stop if all attributes are parsed
        if feature_id and parent_id and biotype:
            break
        key, value = attr.split("=")
        value = value.strip()
        if key == "ID":
            feature_id = value
        elif key == "Parent":
            parent_id = sys.intern(value)
        elif key == "biotype" or key == "gene_biotype" or key == "transcript_biotype":
            biotype = sys.intern(value)

    return feature_type, length, parent_id, biotype, feature_id

FEATURE_CODES = {
    "exon": 1,
    "CDS": 2,
}

feat_code_get = FEATURE_CODES.get

def process_feature(
    feature_id: str,
    feature_type: str,
    length: int,
    parent_id: Optional[str],
    biotype: Optional[str],
    roots: dict,
    id_to_root: dict,
    transcripts: dict,
):
    """Optimized with reduced redundant checks"""
    # Gene (root) feature
    if feature_type in COMMON_GENE_TYPES or (not parent_id and feature_id):
        gene = Gene(feature_type, biotype, length)
        roots[feature_id] = gene
        id_to_root[feature_id] = feature_id
        return True


    processed = False #this is what we return
    
    #getters for cached access
    get_transcript = transcripts.get
    get_root_id = id_to_root.get
    get_root = roots.get
    root_id = get_root_id(parent_id)

    # if root is not found, is an orphan, resolve later
    if root_id is None: 
        return False

    #store feature to root mapping
    if feature_id:
        id_to_root[feature_id] = root_id
    
    feature_code = feat_code_get(feature_type,0) #0 if not an exon or CDS
    # Parent Case 1: The parent is a gene, then this should be a transcript
    if parent_id == root_id:
        if feature_id and (feature_type in COMMON_TRANSCRIPT_TYPES or feature_code == 0):
            t = get_transcript(feature_id)
            if t is None:
                t = Transcript(root_id, feature_type)
                t.length = length
                transcripts[feature_id] = t
            else:
                t.gene_id = root_id
                t.type = feature_type
                t.length = length
            # Mark gene as having children
            get_root(root_id).has_children = True

        processed = True
        return True # exit early if parent is a gene

    # Parent Case 2: The parent is a transcript
    t = get_transcript(parent_id) 
    if t is None:
        t = Transcript(root_id, None)
        transcripts[parent_id] = t
    
    # Mark gene as having children (transcript exists)
    root = get_root(root_id)
    root.has_children = True
    
    # Exon Case
    if feature_code == 1:
        t.exon_count += 1
        #set has_multiple_exons if exon count is 2
        if t.exon_count == 2:
            root.has_multiple_exons = True

        t.exons_lengths.append(length)
        t.exon_len_sum += length
        root.has_exon = True

    # CDS Case
    elif feature_code == 2:
        t.cds_count += 1
        t.cds_len_sum += length
        t.cds_lengths.append(length)
        root.has_cds = True

    processed = True

    return processed


def resolve_orphans(orphans: list[OrphanFeature], roots: dict, id_to_root: dict, transcripts: dict):
    """Optimized orphan resolution"""
    still_orphaned = []

    for orphan in orphans:
        processed = process_feature(
            orphan.feature_id,
            orphan.feature_type,
            orphan.length,
            orphan.parent_id,
            orphan.biotype,
            roots,
            id_to_root,
            transcripts,
        )

        if not processed:
            still_orphaned.append(orphan)

    return still_orphaned


def compute_gene_lengths(roots: dict) -> dict:
    """
    Compute gene lengths for each category.
    Returns a dictionary with categories as keys and arrays of gene lengths as values.
    The arrays are used for median calculation.
    """

    cat_gene_lengths = {
        "coding": array("q"),
        "long_non_coding": array("q"),
        "short_non_coding": array("q"),
        "pseudogene": array("q"),
    }
    for gene_info in roots.values():
        cat = gene_info.category
        if cat:
            cat_gene_lengths[cat].append(gene_info.length) #for length summary (we store a list for median calculation)
    return cat_gene_lengths

def _array_factory(t: str = "q") -> array:
    return array(t)


def init_category_statistics_dictionary() -> dict:
    return {
        "gene_length_list": _array_factory(),
        "types": defaultdict(init_transcript_type_statistics_dictionary),
    }

def init_transcript_type_statistics_dictionary() -> dict:
    return {
        "gene_ids": set(),
        "length_list":  _array_factory(),
        "exon_length_list": _array_factory("i"),
        "spliced_length_list": _array_factory(),
        "intron_length_list": _array_factory(),
        "cds_length_list": _array_factory("i"),
        "protein_length_list": _array_factory(), }


# ============================================================================
# File Handling Utilities
# ============================================================================

def is_url(source: str) -> bool:
    """Check if source is a URL."""
    return isinstance(source, str) and source.startswith(("http://", "https://"))


def is_gzipped(source: str) -> bool:
    """Check if source appears to be gzipped based on extension."""
    return source.endswith(".gz") or source.endswith(".gzip") or source.endswith(".bgz")


def download_url_to_temp(url: str, is_gzipped: bool) -> str:
    """Download URL to temporary file and return the path."""
    fd, temp_path = tempfile.mkstemp(suffix=".gff" if not is_gzipped else ".gff.gz")
    os.close(fd)
    print(f"[gffy] Downloading URL to temp file: {temp_path}")
    with urlopen(url) as r, open(temp_path, "wb") as f:
        shutil.copyfileobj(r, f)
    print(f"[gffy] Downloaded {os.path.getsize(temp_path)} bytes")
    return temp_path


def open_gff_stream(path: str, is_gzipped: bool) -> IO:
    """Open GFF file (gzipped or plain text) as text stream."""
    if is_gzipped:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def stream_gff_lines(file_handle: IO):
    """Stream non-comment, non-empty lines from GFF file."""
    for line in file_handle:
        # Fast path: check first char only (avoid strip() overhead)
        if not line or line[0] == "#" or line[0] in "\n\r":
            continue
        yield line


# ============================================================================
# Attribute Parsing
# ============================================================================

def parse_gff_attributes_fast(attrs: str, feature_code: int) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Fast GFF3 attribute parser with minimal string operations.
    
    Args:
        attrs: The attributes column (9th column) from GFF3 line
        feature_code: 0=transcript, 1=exon, 2=CDS (from FEATURE_CODES)
    
    Returns:
        Tuple of (feature_id, parent_id, biotype)
    """
    intern = sys.intern
    parent_id = None
    biotype = None
    feature_id = None
    
    # Parse attributes without creating list - iterate over string directly
    attr_start = 0
    attrs_len = len(attrs)
    
    while attr_start < attrs_len:
        # Find next semicolon
        semi_pos = attrs.find(";", attr_start)
        if semi_pos == -1:
            semi_pos = attrs_len
        
        # Find '=' in this attribute segment
        eq_pos = attrs.find("=", attr_start, semi_pos)
        if eq_pos != -1:
            key_start = attr_start
            key_end = eq_pos
            val_start = eq_pos + 1
            val_end = semi_pos
            
            # Compare key without extracting substring when possible
            key_len = key_end - key_start
            
            if key_len == 2 and attrs[key_start:key_end] == "ID":
                feature_id = attrs[val_start:val_end].strip()
            elif key_len == 6 and attrs[key_start:key_end] == "Parent":
                parent_id = intern(attrs[val_start:val_end].strip())
                # For exon/CDS, parent is the only required attribute
                if feature_code != 0:  # exon or CDS
                    break
            elif parent_id and (
                (key_len == 7 and attrs[key_start:key_end] == "biotype") or
                (key_len == 12 and attrs[key_start:key_end] == "gene_biotype") or
                (key_len == 17 and attrs[key_start:key_end] == "transcript_biotype")
            ):
                biotype = intern(attrs[val_start:val_end].strip())
                # If we have all needed attributes, stop parsing
                if feature_id and parent_id:
                    break
        
        attr_start = semi_pos + 1
    
    return feature_id, parent_id, biotype


# ============================================================================
# Transcript Processing
# ============================================================================

def process_transcript_feature(
    transcript_id: str,
    parent_id: str,
    feature_type: str,
    biotype: Optional[str],
    length: int,
    transcripts: dict
) -> None:
    """Process a transcript feature and update transcripts dictionary."""
    transcript = transcripts.get(transcript_id)
    if transcript is None:
        transcript = Transcript(parent_id, feature_type)
        transcripts[transcript_id] = transcript
    else:
        # Update existing transcript
        transcript.gene_id = sys.intern(parent_id)
        transcript.type = sys.intern(feature_type) if feature_type else None
    
    transcript.length = length


def process_exon_or_cds(
    parent_id: str,
    feature_code: int,
    length: int,
    transcripts: dict
) -> None:
    """Process exon or CDS feature and update parent transcript."""
    transcript = transcripts.get(parent_id)
    if transcript is None:
        transcript = Transcript(parent_id, None)
        transcripts[parent_id] = transcript
    
    if feature_code == 1:  # exon
        transcript.exon_count += 1
        transcript.exons_lengths.append(length)
        transcript.exon_len_sum += length
    else:  # CDS (feature_code == 2)
        transcript.cds_count += 1
        transcript.cds_lengths.append(length)
        transcript.cds_len_sum += length
