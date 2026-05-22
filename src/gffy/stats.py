"""
Fast, low-RAM GFF3 statistics calculator.

Output schema matches annotrieve ``compute_features_statistics`` (GFFStats).
Small files use a single in-memory pass with slim records; large files are
sorted by seqid into a temporary file and processed per-seqid to cap peak RAM.
"""

from __future__ import annotations

import gzip
import os
import sys
import tempfile
from array import array
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Iterable, Iterator
from ._io import is_gzipped_path, is_url, iter_text_lines
from .cache import (
    urlopen_request,
    SourceInfo,
    build_sorted_gff,
    inspect_source,
    is_big,
)

MISSING_BIOTYPE = "biotype_missing"

GENE_CODES = frozenset({
    "gene",
    "ncRNA_gene",
    "pseudogene",
})

TRANSCRIPT_CODES = frozenset({
    "mRNA",
    "ncRNA",
    "tRNA",
    "rRNA",
    "sRNA",
    "snRNA",
    "snoRNA",
    "misc_RNA",
    "miRNA",
    "piRNA",
    "siRNA",
    "lncRNA",
    "transcript",
})

SUB_FEATURE_CODES = frozenset({
    "exon",
    "CDS",
})

DNA_REGION_CODES = frozenset({
    "chromosome",
    "contig",
    "scaffold",
    "plasmid",
    "mitochondrion",
    "chloroplast",
    "mitogenome",
    "region",
    "biological_region",
    "telomere",
    "supercontig",
    "scaffold_region",
    "repeat_region",
    "mobile_genetic_element"
})


def parse_attributes_full(attr_string: str) -> dict[str, str]:
    """Parse all GFF3 attributes from column 9."""
    d: dict[str, str] = {}
    if not attr_string:
        return d
    for part in attr_string.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = v.strip()
    return d


class _FeatureLengthStats:
    __slots__ = ("total_count", "mean_length", "min_length", "max_length", "sum_lengths")

    def __init__(self) -> None:
        self.total_count = 0
        self.mean_length = 0.0
        self.min_length: int | None = None
        self.max_length: int | None = None
        self.sum_lengths = 0

    def update_length(self, length: int) -> None:
        self.total_count += 1
        self.sum_lengths += length
        if self.min_length is None or length < self.min_length:
            self.min_length = length
        if self.max_length is None or length > self.max_length:
            self.max_length = length
        self.mean_length = self.sum_lengths / self.total_count


class _GeneCategoryStats:
    __slots__ = ("count", "length_stats", "biotype_counts", "transcript_type_counts")

    def __init__(self) -> None:
        self.count = 0
        self.length_stats = _FeatureLengthStats()
        self.biotype_counts: dict[str, int] = defaultdict(int)
        self.transcript_type_counts: dict[str, int] = defaultdict(int)


class _TranscriptStats:
    __slots__ = (
        "count",
        "transcript_lengths",
        "exon_counts",
        "cds_counts",
        "concat_exon_lengths",
        "concat_cds_lengths",
        "genes_with_this_type",
        "has_multiple_exons",
        "has_cds",
        "biotype_counts",
        "gene_categories",
    )

    def __init__(self) -> None:
        self.count = 0
        self.transcript_lengths = _FeatureLengthStats()
        self.exon_counts = _FeatureLengthStats()
        self.cds_counts = _FeatureLengthStats()
        self.concat_exon_lengths = _FeatureLengthStats()
        self.concat_cds_lengths = _FeatureLengthStats()
        self.genes_with_this_type: set[str] = set()
        self.has_multiple_exons = False
        self.has_cds = False
        self.biotype_counts: dict[str, int] = defaultdict(int)
        self.gene_categories: dict[str, set[str]] = defaultdict(set)


class _GeneRecord:
    __slots__ = ("feature_type", "biotype", "length", "category")

    def __init__(self, feature_type: str, biotype: str, length: int) -> None:
        self.feature_type = feature_type
        self.biotype = biotype
        self.length = length
        self.category: str | None = None


class _SlimFeature:
    """Transcript-like or alternate-gene feature without storing full attr dict."""

    __slots__ = (
        "feature_type",
        "length",
        "parent_id",
        "transcript_biotype",
        "gene_biotype",
    )

    def __init__(
        self,
        feature_type: str,
        length: int,
        parent_id: str | None,
        transcript_biotype: str | None,
        gene_biotype: str | None,
    ) -> None:
        self.feature_type = feature_type
        self.length = length
        self.parent_id = parent_id
        self.transcript_biotype = transcript_biotype
        self.gene_biotype = gene_biotype


class _GlobalState:
    __slots__ = ("gene_categories", "transcript_stats", "known_transcript_types")

    def __init__(self) -> None:
        self.gene_categories = {
            "coding": _GeneCategoryStats(),
            "pseudogene": _GeneCategoryStats(),
            "non_coding": _GeneCategoryStats(),
        }
        self.transcript_stats: dict[str, _TranscriptStats] = defaultdict(_TranscriptStats)
        self.known_transcript_types: set[str] = set()


def _length_stats_dict(stats: _FeatureLengthStats) -> dict:
    return {
        "min": stats.min_length if stats.min_length is not None else 0,
        "max": stats.max_length if stats.max_length is not None else 0,
        "mean": round(stats.mean_length if stats.mean_length > 0 else 0.0, 2),
    }


def _intern_optional(value: str | None) -> str | None:
    return sys.intern(value) if value else None


def _parse_slim_fields(attr: dict[str, str]) -> tuple[str | None, str | None, str | None]:
    parent = attr.get("Parent") or attr.get("gene") or attr.get("Gene")
    transcript_biotype = attr.get("biotype") or attr.get("transcript_biotype")
    gene_biotype = attr.get("biotype") or attr.get("gene_biotype") or attr.get("type")
    return (
        _intern_optional(parent.strip()) if parent else None,
        _intern_optional(transcript_biotype) if transcript_biotype else None,
        _intern_optional(gene_biotype.lower()) if gene_biotype else None,
    )


def _classify_line(
    line: str,
    exons_cds: dict,
    genes_by_id: dict[str, _GeneRecord],
    features_by_id: dict[str, _SlimFeature],
) -> None:
    """Route one GFF line into per-batch accumulators (slim records only)."""
    f = line.split("\t", 8)
    if len(f) < 9:
        return

    feature_type = sys.intern(f[2])

    if feature_type in SUB_FEATURE_CODES:
        try:
            start, end = int(f[3]), int(f[4])
            length = end - start + 1
        except (ValueError, IndexError):
            return

        attr = parse_attributes_full(f[8])
        parent_ids = attr.get("Parent", "")
        if not parent_ids:
            return

        for parent_id in parent_ids.split(","):
            parent_id = parent_id.strip()
            if not parent_id:
                continue
            parent_id = sys.intern(parent_id)
            if parent_id not in exons_cds:
                exons_cds[parent_id] = {
                    "exon_lengths": array("i"),
                    "cds_lengths": array("i"),
                }
            if feature_type == "exon":
                exons_cds[parent_id]["exon_lengths"].append(length)
            else:
                exons_cds[parent_id]["cds_lengths"].append(length)
        return

    if feature_type in DNA_REGION_CODES:
        return

    try:
        start, end = int(f[3]), int(f[4])
        length = end - start + 1
    except (ValueError, IndexError):
        return

    attr = parse_attributes_full(f[8])
    feature_id = attr.get("ID")
    if not feature_id:
        return

    if feature_type in GENE_CODES:
        gene_ftype = sys.intern(feature_type.lower())
        raw_biotype = attr.get("biotype") or attr.get("gene_biotype") or attr.get("type")
        gene_biotype = (
            sys.intern(raw_biotype.lower()) if raw_biotype else MISSING_BIOTYPE
        )
        if feature_id not in genes_by_id:
            genes_by_id[feature_id] = _GeneRecord(gene_ftype, gene_biotype, length)
        return

    parent_id, transcript_biotype, gene_biotype = _parse_slim_fields(attr)
    if feature_id not in features_by_id:
        features_by_id[feature_id] = _SlimFeature(
            feature_type, length, parent_id, transcript_biotype, gene_biotype
        )


def _resolve_batch(
    exons_cds: dict,
    genes_by_id: dict[str, _GeneRecord],
    features_by_id: dict[str, _SlimFeature],
    state: _GlobalState,
) -> None:
    """Drain one batch (whole file or single seqid) into global aggregates."""
    gene_categories = state.gene_categories
    transcript_stats = state.transcript_stats
    known_transcript_types = state.known_transcript_types

    transcripts: dict = {}

    for tid, rec in features_by_id.items():
        if tid not in exons_cds:
            continue
        if rec.feature_type in DNA_REGION_CODES or rec.feature_type in GENE_CODES:
            continue
        if rec.feature_type in SUB_FEATURE_CODES:
            continue

        exon_lengths = exons_cds[tid]["exon_lengths"]
        cds_lengths = exons_cds[tid]["cds_lengths"]
        if len(exon_lengths) == 0 and len(cds_lengths) == 0:
            continue

        # Match annotrieve feature_stats: keep non-DNA feature types (e.g. lnc_RNA).
        if rec.feature_type in TRANSCRIPT_CODES:
            ts_type = rec.feature_type
        elif rec.feature_type not in DNA_REGION_CODES:
            ts_type = rec.feature_type
        else:
            ts_type = sys.intern("transcript")

        transcripts[tid] = {
            "type": ts_type,
            "biotype": rec.transcript_biotype,
            "gene": rec.parent_id,
            "length": rec.length,
            "exon_lengths": exon_lengths,
            "cds_lengths": cds_lengths,
        }

    gene_has_exon: set[str] = set()
    gene_has_cds: set[str] = set()

    for tid, tdata in transcripts.items():
        gene_id = tdata.get("gene")
        if gene_id:
            if len(tdata.get("exon_lengths", array("i"))) > 0:
                gene_has_exon.add(gene_id)
            if len(tdata.get("cds_lengths", array("i"))) > 0:
                gene_has_cds.add(gene_id)

    for fid, rec in features_by_id.items():
        if rec.feature_type in TRANSCRIPT_CODES or rec.feature_type in SUB_FEATURE_CODES:
            continue
        if rec.feature_type in DNA_REGION_CODES or rec.feature_type in GENE_CODES:
            continue
        if fid not in gene_has_exon and fid not in gene_has_cds:
            continue
        gene_ftype = sys.intern(rec.feature_type.lower())
        gene_biotype = rec.gene_biotype or MISSING_BIOTYPE
        if fid not in genes_by_id:
            genes_by_id[fid] = _GeneRecord(gene_ftype, gene_biotype, rec.length)

    for gene_id, info in genes_by_id.items():
        ftype = info.feature_type
        biotype = info.biotype
        length = info.length
        has_cds = gene_id in gene_has_cds
        has_exon = gene_id in gene_has_exon

        if ftype == "pseudogene":
            category = "pseudogene"
        elif has_cds or biotype == "protein_coding":
            category = "coding"
        elif has_exon:
            category = "non_coding"
        else:
            continue

        gene_categories[category].count += 1
        gene_categories[category].length_stats.update_length(length)
        biotype_key = biotype if biotype else MISSING_BIOTYPE
        gene_categories[category].biotype_counts[biotype_key] += 1
        info.category = category

    for tid, tdata in transcripts.items():
        ts_type = tdata.get("type", "transcript")
        ts_biotype = tdata.get("biotype", "")
        ts_gene = tdata.get("gene")
        known_transcript_types.add(ts_type)
        ts = transcript_stats[ts_type]

        exon_lengths = tdata.get("exon_lengths", array("i"))
        cds_lengths = tdata.get("cds_lengths", array("i"))
        transcript_length = tdata.get("length", 0)
        exon_count = len(exon_lengths)
        cds_count = len(cds_lengths)

        if exon_count > 1:
            ts.has_multiple_exons = True
        if cds_count > 0:
            ts.has_cds = True

        for exon_len in exon_lengths:
            ts.exon_counts.update_length(exon_len)
        for cds_len in cds_lengths:
            ts.cds_counts.update_length(cds_len)

        concat_exon_len = sum(exon_lengths)
        if concat_exon_len > 0:
            ts.concat_exon_lengths.update_length(concat_exon_len)

        concat_cds_len = sum(cds_lengths)
        if concat_cds_len > 0:
            ts.concat_cds_lengths.update_length(concat_cds_len)

        ts.transcript_lengths.update_length(transcript_length)
        ts.count += 1

        if ts_gene:
            ts.genes_with_this_type.add(ts_gene)
            if ts_gene in genes_by_id:
                gene_category = genes_by_id[ts_gene].category
                if gene_category:
                    ts.gene_categories[gene_category].add(ts_gene)
                    gene_categories[gene_category].transcript_type_counts[ts_type] += 1

        biotype_key = ts_biotype if ts_biotype else MISSING_BIOTYPE
        ts.biotype_counts[biotype_key] += 1


def _finalize(state: _GlobalState) -> dict:
    gene_category_stats_dict: dict = {}
    for category in ("coding", "pseudogene", "non_coding"):
        stats = state.gene_categories[category]
        if stats.count > 0:
            gene_category_stats_dict[category] = {
                "total_count": stats.count,
                "length_stats": _length_stats_dict(stats.length_stats),
                "biotype_counts": dict(sorted(stats.biotype_counts.items())),
                "transcript_type_counts": dict(
                    sorted(stats.transcript_type_counts.items())
                ),
            }

    transcript_type_stats_dict: dict = {}
    sorted_transcript_types = sorted(
        [
            ttype
            for ttype in state.known_transcript_types
            if ttype in state.transcript_stats
            and state.transcript_stats[ttype].count > 0
        ],
        key=lambda ttype: state.transcript_stats[ttype].count,
        reverse=True,
    )

    for ttype in sorted_transcript_types:
        ts = state.transcript_stats[ttype]
        tl = ts.transcript_lengths
        ec = ts.exon_counts
        cel = ts.concat_exon_lengths

        exon_stats_data: dict = {
            "length": _length_stats_dict(ec),
            "total_count": ec.total_count,
        }
        if ec.total_count > ts.count:
            exon_stats_data["concatenated_length"] = _length_stats_dict(cel)

        entry: dict = {
            "length_stats": _length_stats_dict(tl),
            "total_count": ts.count,
            "biotype_counts": dict(sorted(ts.biotype_counts.items())),
            "associated_genes": {
                "total_count": len(ts.genes_with_this_type),
                "gene_categories": dict(
                    sorted((cat, len(gene_set)) for cat, gene_set in ts.gene_categories.items())
                ),
            },
            "exon_stats": exon_stats_data,
        }

        if ts.has_cds:
            cc = ts.cds_counts
            ccdl = ts.concat_cds_lengths
            entry["cds_stats"] = {
                "total_count": cc.total_count,
                "length": _length_stats_dict(cc),
                "concatenated_length": _length_stats_dict(ccdl),
            }

        transcript_type_stats_dict[ttype] = entry

    return {
        "gene_category_stats": gene_category_stats_dict,
        "transcript_type_stats": transcript_type_stats_dict,
    }


def _empty_scratch() -> tuple[dict, dict, dict]:
    return {}, {}, {}


def _iter_gff_lines(source: str) -> Iterator[str]:
    """Yield non-comment GFF lines from a local path or URL (plain or .gz)."""
    if is_url(source):
        raw = urlopen_request(source, timeout=300)
        if is_gzipped_path(source):
            with gzip.GzipFile(fileobj=raw) as gz:
                for raw_line in gz:
                    line = raw_line.decode("utf-8", errors="replace")
                    if line.startswith("#") or not line.strip():
                        continue
                    yield line.rstrip("\n")
        else:
            for raw_line in raw:
                line = raw_line.decode("utf-8", errors="replace")
                if line.startswith("#") or not line.strip():
                    continue
                yield line.rstrip("\n")
        return

    yield from iter_text_lines(source)


def _compute_per_seqid(sorted_path: str | Path) -> dict:
    state = _GlobalState()
    current_seqid: str | None = None
    exons_cds, genes_by_id, features_by_id = _empty_scratch()

    for line in iter_text_lines(sorted_path):
        tab = line.find("\t")
        if tab <= 0:
            continue
        seqid = line[:tab]

        if current_seqid is not None and seqid != current_seqid:
            _resolve_batch(exons_cds, genes_by_id, features_by_id, state)
            exons_cds, genes_by_id, features_by_id = _empty_scratch()

        current_seqid = seqid
        _classify_line(line, exons_cds, genes_by_id, features_by_id)

    if current_seqid is not None:
        _resolve_batch(exons_cds, genes_by_id, features_by_id, state)

    return _finalize(state)


def compute_gff_stats_from_lines(lines: Iterable[str]) -> dict:
    """Compute statistics from an iterable of GFF lines (unordered-safe)."""
    exons_cds, genes_by_id, features_by_id = _empty_scratch()
    for line in lines:
        _classify_line(line, exons_cds, genes_by_id, features_by_id)
    state = _GlobalState()
    _resolve_batch(exons_cds, genes_by_id, features_by_id, state)
    return _finalize(state)


def _compute_low_memory(source: str, info: SourceInfo) -> dict:
    """Sort into a temp file, compute per-seqid, then discard the temp dir."""
    with tempfile.TemporaryDirectory(prefix="gffy-") as tmp:
        sorted_out = Path(tmp) / "sorted.gff.gz"
        build_sorted_gff(source, sorted_out, info)
        return _compute_per_seqid(sorted_out)


def compute_gff_stats(
    source: str,
    *,
    force_low_memory: bool = False,
) -> dict:
    """
    Compute GFF3 feature statistics from a local path or URL.

    Large files (gz > 100 MB or plain > 1 GB) or ``force_low_memory=True`` sort
    into a temporary file, stream per-seqid, then remove the temp file.
    """
    info = inspect_source(source)
    if force_low_memory or is_big(info):
        return _compute_low_memory(source, info)
    return compute_gff_stats_from_lines(_iter_gff_lines(source))


def describe_compute_mode(
    source: str,
    *,
    force_low_memory: bool = False,
) -> tuple[str, SourceInfo]:
    """Return human-readable mode label and source metadata for logging."""
    info = inspect_source(source)
    if force_low_memory or is_big(info):
        return "low-memory (per-seqid)", info
    return "single-pass", info
