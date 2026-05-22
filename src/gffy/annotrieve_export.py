"""
Build Annotrieve Favorites-compatible CustomAnnotation records from GFF3.

Uses sorted GFF (matching server upload pipeline) for MD5 identity, then two passes:
pass 1 — MD5 + features_summary in one stream; pass 2 — features_statistics.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import tempfile
from urllib.parse import urlparse
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

from .cache import build_sorted_gff, inspect_source, is_big
from ._io import is_gzipped_path, is_url, iter_text_lines
from .stats import _compute_per_seqid, compute_gff_stats_from_lines

def derive_custom_name(source: str, override: str | None = None) -> str:
    """
    Derive a display name from a path or URL.

    For http(s)/ftp URLs, uses only the URL path (host/domain omitted).
    Replaces ``/`` and spaces with ``_``, collapses repeats.
    *override* wins when non-empty (single-file CLI).
    """
    if override and override.strip():
        return override.strip()
    raw = source.strip()
    if is_url(raw):
        name = urlparse(raw).path
    else:
        name = raw
    name = name.replace("/", "_").replace(" ", "_")
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "custom_annotation"


class _SummaryState:
    __slots__ = (
        "attribute_keys",
        "types",
        "sources",
        "biotypes",
        "root_type_counts",
        "types_missing_id",
    )

    def __init__(self) -> None:
        self.attribute_keys: set[str] = set()
        self.types: set[str] = set()
        self.sources: set[str] = set()
        self.biotypes: set[str] = set()
        self.root_type_counts: dict[str, int] = {}
        self.types_missing_id: set[str] = set()

    def finalize(self) -> dict[str, Any]:
        has_cds = "CDS" in self.types
        has_exon = "exon" in self.types
        has_biotype = bool(self.biotypes)
        return {
            "attribute_keys": sorted(self.attribute_keys),
            "types": sorted(self.types),
            "sources": sorted(self.sources),
            "biotypes": sorted(self.biotypes),
            "root_type_counts": dict(self.root_type_counts),
            "types_missing_id": sorted(self.types_missing_id),
            "has_biotype": has_biotype,
            "has_cds": has_cds,
            "has_exon": has_exon,
        }


def _accumulate_summary_from_line(line: str, state: _SummaryState) -> None:
    """Port of annotrieve feature_summary per-line logic."""
    fields = line.split("\t")
    if len(fields) < 9:
        return

    source = fields[1]
    feature_type = fields[2]
    attributes = fields[8]

    has_id = False
    has_parent = False

    for part in attributes.split(";"):
        part = part.strip()
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key == "ID":
            has_id = True
        elif key == "Parent":
            has_parent = True
        state.attribute_keys.add(key)
        if key in ("biotype", "gene_biotype", "transcript_biotype"):
            state.biotypes.add(value)

    if not has_id:
        state.types_missing_id.add(feature_type)
    if not has_parent:
        state.root_type_counts[feature_type] = (
            state.root_type_counts.get(feature_type, 0) + 1
        )

    state.types.add(feature_type)
    state.sources.add(source)


def compute_features_summary_from_lines(lines: Iterator[str]) -> dict[str, Any]:
    """Compute features_summary from an iterable of GFF data lines (no comments)."""
    state = _SummaryState()
    for line in lines:
        _accumulate_summary_from_line(line, state)
    return state.finalize()


def hash_and_features_summary_from_sorted_gz(path: str | Path) -> tuple[str, int, dict[str, Any]]:
    """
    Pass 1: MD5 over decompressed sorted bytes + features_summary in one read.

    *path* may be plain ``.gff`` or ``.gff.gz`` (sorted seqid order).
    """
    path = Path(path)
    hasher = hashlib.md5()
    uncompressed_size = 0
    state = _SummaryState()

    if is_gzipped_path(path):
        opener = gzip.open
        mode = "rb"
    else:
        opener = open  # type: ignore[assignment]
        mode = "rb"

    with opener(path, mode) as fh:
        for raw_line in fh:
            hasher.update(raw_line)
            uncompressed_size += len(raw_line)
            if raw_line.startswith(b"#") or not raw_line.strip():
                continue
            line = raw_line.decode("utf-8", errors="replace").rstrip("\n\r")
            _accumulate_summary_from_line(line, state)

    return hasher.hexdigest(), uncompressed_size, state.finalize()


@contextmanager
def prepare_sorted_gff_gz(source: str):
    """
    Yield path to a sorted gzipped GFF in a temporary directory.

    Always sorts for Annotrieve export (MD5 must match sorted content).
    """
    info = inspect_source(source)
    with tempfile.TemporaryDirectory(prefix="gffy-annotrieve-") as tmp:
        sorted_gz = Path(tmp) / "sorted.gff.gz"
        build_sorted_gff(source, sorted_gz, info)
        yield sorted_gz


def build_custom_annotation(
    source: str,
    *,
    custom_name_override: str | None = None,
    force_low_memory: bool = False,
) -> dict[str, Any]:
    """
    Build a CustomAnnotation-compatible dict for Annotrieve Favorites import.

    Sorts the GFF, then pass 1 (hash + summary) and pass 2 (statistics).

    Pass 2 follows ``compute_gff_stats`` policy: per-seqid only when the source
    is large or ``force_low_memory``; otherwise single-pass over the sorted file.
    """
    custom_name = derive_custom_name(source, custom_name_override)
    info = inspect_source(source)
    use_per_seqid = force_low_memory or is_big(info)

    with prepare_sorted_gff_gz(source) as sorted_gz:
        md5_hex, file_size, features_summary = hash_and_features_summary_from_sorted_gz(
            sorted_gz
        )
        if not features_summary.get("types") or not features_summary.get("sources"):
            raise ValueError("Annotation has no types or sources, skipping...")

        if use_per_seqid:
            features_statistics = _compute_per_seqid(sorted_gz)
        else:
            features_statistics = compute_gff_stats_from_lines(
                iter_text_lines(sorted_gz)
            )

    uploaded_at = datetime.now(timezone.utc).isoformat()

    return {
        "kind": "custom",
        "annotation_id": md5_hex,
        "custom_name": custom_name,
        "uploaded_md5": md5_hex,
        "uploaded_at": uploaded_at,
        "uploaded_file_size": file_size,
        "features_summary": features_summary,
        "features_statistics": features_statistics,
    }


def write_annotrieve_json(record: dict[str, Any], path: str | Path) -> None:
    """Write one CustomAnnotation JSON file (pretty-printed)."""
    path = Path(path)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2)
        fh.write("\n")


def format_annotrieve_jsonl_line(record: dict[str, Any]) -> str:
    """Serialize one record as a compact JSONL line."""
    return json.dumps(record, separators=(",", ":"))
