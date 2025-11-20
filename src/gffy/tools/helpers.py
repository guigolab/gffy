from collections import defaultdict
from typing import Optional
from .classes import Gene, Transcript, OrphanFeature
import sys
import requests
import gzip
from contextlib import contextmanager
from array import array

@contextmanager
def open_stream(path, chunk_size=1024, is_gzipped=True):
    """
    Yield a line-by-line stream for local or HTTP/HTTPS files.
    Supports gzip automatically.
    """
    is_url = path.startswith(("http://", "https://"))

    if is_url:
        resp = requests.get(path, stream=True)
        resp.raise_for_status()

        # Check gzip by header OR file extension
        http_is_gzip = (
            resp.headers.get("Content-Encoding") == "gzip"
            or is_gzipped
        )

        if http_is_gzip:
            # Wrap raw response in a GzipFile for streaming decompression
            gz = gzip.GzipFile(fileobj=resp.raw)

            def line_stream():
                for line in gz:
                    yield line.decode("utf-8")

            try:
                yield line_stream()
            finally:
                gz.close()
                resp.close()

        else:
            # Plain text streaming
            def line_stream():
                for line in resp.iter_lines(chunk_size=chunk_size, decode_unicode=True):
                    if line:
                        yield line

            try:
                yield line_stream()
            finally:
                resp.close()

    else:
        # Local file path
        if is_gzipped:
            f = gzip.open(path, "rt", encoding="utf-8")
        else:
            f = open(path, "r", encoding="utf-8")

        try:
            yield f
        finally:
            f.close()

def init_orphan_feature(
    feature_id: str,
    feature_type: str,
    length: int,
    parent_ids: list[str],
    biotype: Optional[str],
) -> OrphanFeature:
    return OrphanFeature(feature_id, feature_type, length, parent_ids, biotype)

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
    parent_ids = []
    biotype = None
    feature_id = None
    for attr in attr_str.split(";"):
        key, value = attr.split("=")
        if key == "ID":
            feature_id = value.strip()
        elif key == "Parent":
            parent_ids = [p.strip() for p in (value.split(",") if "," in value else [value])]
        elif key == "biotype" or key == "gene_biotype" or key == "transcript_biotype":
            biotype = sys.intern(value.strip())

    return feature_type, length, parent_ids, biotype, feature_id

FEATURE_CODES = {
    "exon": 1,
    "CDS": 2,
}

feat_code_get = FEATURE_CODES.get

def process_feature(
    feature_id: str,
    feature_type: str,
    length: int,
    parent_ids: list[str],
    biotype: Optional[str],
    roots: dict,
    id_to_root: dict,
    transcripts: dict,
):
    """Optimized with reduced redundant checks"""

    # Gene feature
    if not parent_ids and feature_id:
        gene = Gene(feature_type, biotype, length)
        roots[feature_id] = gene
        id_to_root[feature_id] = feature_id
        return True

    processed = False #this is what we return

    feature_code = feat_code_get(feature_type,0) #0 if not an exon or CDS

    #getters for quick access
    get_transcript = transcripts.get
    get_root_id = id_to_root.get
    get_root = roots.get
    for parent_id in parent_ids:
        root_id = get_root_id(parent_id)

        # if root is not found, is an orphan, resolve later
        if root_id is None: 
            continue

        #store feature to root mapping
        if feature_id:
            id_to_root[feature_id] = root_id
        
        # Parent Case 1: The parent is a gene, then this is a transcript
        if parent_id == root_id:
            if feature_id and feature_code == 0:
                t = get_transcript(feature_id) #transcripts.get(feature_id)
                if t is None:
                    t = Transcript(root_id, feature_type)
                    t.length = length
                    transcripts[feature_id] = t
                else:
                    t.gene_id = root_id
                    t.type = feature_type
                    t.length = length

            processed = True
            continue # exit early if parent is a gene

        # Parent Case 2: The parent is a transcript
        t = get_transcript(parent_id) 
        if t is None:
            t = Transcript(root_id, None)
            transcripts[parent_id] = t
        
        # Exon Case
        if feature_code == 1:
            t.exon_count += 1
            #set has_multiple_exons if exon count is 2
            if t.exon_count == 2:
                get_root(root_id).has_multiple_exons = True

            t.exons_lengths.append(length)
            t.exon_len_sum += length
            get_root(root_id).has_exon = True

        # CDS Case
        elif feature_code == 2:
            t.cds_count += 1
            t.cds_len_sum += length
            t.cds_lengths.append(length)
            get_root(root_id).has_cds = True

        processed = True

    return processed


def resolve_orphans(orphans: list[OrphanFeature], roots: dict, id_to_root: dict, transcripts: dict):
    """Optimized orphan resolution"""
    still_orphaned = []

    for orphan in orphans:
        length = orphan.end - orphan.start + 1
        processed = process_feature(
            orphan.feature_id,
            orphan.feature_type,
            length,
            orphan.parent_ids,
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
