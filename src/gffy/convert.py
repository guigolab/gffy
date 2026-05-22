"""
Fast, low-RAM GTF to GFF3 converter.

Translates GTF attribute nomenclature (gene_id, transcript_id) to GFF3 ID/Parent,
infers biotype=protein_coding on root features with CDS-bearing transcripts, and
rewrites transcript rows to mRNA when they have CDS children.
"""

from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path
from typing import Iterator, TextIO

from ._io import is_gzipped_path, iter_text_lines, open_text_auto
from .cache import (
    SourceInfo,
    _download_url,
    _sort_to_gzip,
    inspect_source,
    is_big,
)

SUB_FEATURE_TYPES = frozenset({
    "exon",
    "CDS",
    "start_codon",
    "stop_codon",
    "five_prime_utr",
    "five_prime_UTR",
    "three_prime_utr",
    "three_prime_UTR",
    "UTR",
    "Selenocysteine",
    "intron",
})

_GTF_ATTR_RE = re.compile(
    r'([^\s=;]+)\s+(?:"([^"]*)"|\'([^\']*)\'|([^;\s]+))'
)

_GFF3_ENCODE = str.maketrans({
    ";": "%3B",
    "=": "%3D",
    "&": "%26",
    ",": "%2C",
    "\t": "%09",
    "\r": "%0D",
    "\n": "%0A",
    "%": "%25",
})


def parse_gtf_attributes(attr_string: str) -> list[tuple[str, str]]:
    """Parse GTF column 9 into ordered (key, value) pairs."""
    pairs: list[tuple[str, str]] = []
    s = attr_string.strip().rstrip(";")
    if not s:
        return pairs
    for m in _GTF_ATTR_RE.finditer(s):
        key = m.group(1)
        val = m.group(2) or m.group(3) or m.group(4) or ""
        pairs.append((key, val))
    return pairs


def _encode_gff3_value(value: str) -> str:
    return value.translate(_GFF3_ENCODE)


def _format_gff3_attrs(pairs: list[tuple[str, str]]) -> str:
    """Emit GFF3 column 9 with ID first, Parent second."""
    id_pair: tuple[str, str] | None = None
    parent_pair: tuple[str, str] | None = None
    rest: list[tuple[str, str]] = []
    for key, val in pairs:
        if key == "ID":
            id_pair = (key, val)
        elif key == "Parent":
            parent_pair = (key, val)
        else:
            rest.append((key, val))
    ordered: list[tuple[str, str]] = []
    if id_pair:
        ordered.append(id_pair)
    if parent_pair:
        ordered.append(parent_pair)
    ordered.extend(rest)
    return ";".join(f"{k}={_encode_gff3_value(v)}" for k, v in ordered)


def _scan_pass1(
    lines: Iterator[str],
) -> tuple[set[str], set[str]]:
    """
    First pass: find transcripts with CDS and genes that inherit CDS.

    Returns (transcripts_with_cds, genes_with_cds).
    """
    transcripts_with_cds: set[str] = set()
    genes_with_cds: set[str] = set()
    transcript_to_gene: dict[str, str] = {}

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        ftype = parts[2]
        attrs = dict(parse_gtf_attributes(parts[8]))
        gene_id = attrs.get("gene_id")
        transcript_id = attrs.get("transcript_id")

        if transcript_id and gene_id:
            transcript_to_gene[sys.intern(transcript_id)] = sys.intern(gene_id)

        if ftype == "CDS" and transcript_id:
            tid = sys.intern(transcript_id)
            transcripts_with_cds.add(tid)
            gid = transcript_to_gene.get(tid)
            if gid:
                genes_with_cds.add(gid)
            elif gene_id:
                genes_with_cds.add(sys.intern(gene_id))

    for tid in transcripts_with_cds:
        gid = transcript_to_gene.get(tid)
        if gid:
            genes_with_cds.add(gid)

    return transcripts_with_cds, genes_with_cds


def _convert_line(
    line: str,
    transcripts_with_cds: set[str],
    genes_with_cds: set[str],
) -> str | None:
    """Convert one GTF data line to GFF3. Returns None to skip."""
    parts = line.split("\t")
    if len(parts) < 9:
        return None

    seqid, source, ftype, start, end, score, strand, phase, attr_col = parts[:9]
    pairs = parse_gtf_attributes(attr_col)
    attrs = dict(pairs)
    gene_id = attrs.get("gene_id")
    transcript_id = attrs.get("transcript_id")
    exon_id = attrs.get("exon_id")

    out_type = ftype
    if ftype == "transcript" and transcript_id and transcript_id in transcripts_with_cds:
        out_type = "mRNA"

    is_sub = ftype in SUB_FEATURE_TYPES
    is_root = bool(gene_id) and (
        not transcript_id or transcript_id == gene_id
    )
    is_transcript_level = (
        bool(transcript_id)
        and bool(gene_id)
        and transcript_id != gene_id
        and not is_sub
    )

    out_pairs: list[tuple[str, str]] = []
    injected_biotype = False

    if is_root and gene_id:
        out_pairs.append(("ID", gene_id))
        if gene_id in genes_with_cds and "biotype" not in attrs:
            out_pairs.append(("biotype", "protein_coding"))
            injected_biotype = True
    elif is_transcript_level:
        out_pairs.append(("ID", transcript_id))  # type: ignore[arg-type]
        out_pairs.append(("Parent", gene_id))  # type: ignore[arg-type]
    elif is_sub and transcript_id:
        out_pairs.append(("Parent", transcript_id))
        if exon_id and ftype == "exon":
            out_pairs.append(("ID", exon_id))

    skip_keys = frozenset({"gene_id", "transcript_id", "exon_id"})
    for key, val in pairs:
        if key in skip_keys:
            continue
        if key == "biotype" and injected_biotype:
            continue
        out_pairs.append((key, val))

    col9 = _format_gff3_attrs(out_pairs)
    return "\t".join(
        (seqid, source, out_type, start, end, score, strand, phase, col9)
    )


def _emit_pass2(
    lines: Iterator[str],
    out: TextIO,
    transcripts_with_cds: set[str],
    genes_with_cds: set[str],
) -> int:
    """Second pass: write GFF3 lines. Returns feature count."""
    out.write("##gff-version 3\n")
    count = 0
    for line in lines:
        converted = _convert_line(line, transcripts_with_cds, genes_with_cds)
        if converted:
            out.write(converted + "\n")
            count += 1
    return count


def _prepare_input_path(
    source: str,
    info: SourceInfo,
    tmp: Path,
    *,
    force_low_memory: bool,
    sort: bool,
) -> Path:
    """Stage source under *tmp* and optionally sort by seqid; return path to read."""
    if info.is_url:
        suffix = ".gtf.gz" if info.is_gzipped else ".gtf"
        local = tmp / f"raw{suffix}"
        _download_url(source, local)
    else:
        local = Path(info.local_path or source)

    use_sort = sort and (force_low_memory or is_big(info))
    if use_sort:
        sorted_out = tmp / "sorted.gtf.gz"
        is_gz = info.is_gzipped or is_gzipped_path(local)
        _sort_to_gzip(local, sorted_out, is_gz)
        return sorted_out

    return local


def _run_conversion(input_path: Path, output: Path) -> tuple[int, int, int]:
    """Scan and convert *input_path* GTF into *output* GFF3."""
    transcripts_with_cds, genes_with_cds = _scan_pass1(iter_text_lines(input_path))
    with open_text_auto(output, "wt") as out_fh:
        feature_count = _emit_pass2(
            iter_text_lines(input_path),
            out_fh,
            transcripts_with_cds,
            genes_with_cds,
        )
    return feature_count, len(genes_with_cds), len(transcripts_with_cds)


def describe_convert_mode(
    source: str,
    *,
    force_low_memory: bool = False,
    sort: bool = True,
) -> tuple[str, SourceInfo]:
    """Return human-readable mode label and source metadata for logging."""
    info = inspect_source(source)
    if sort and (force_low_memory or is_big(info)):
        return "low-memory (sorted staging)", info
    return "single-pass", info


def convert_gtf_to_gff3(
    source: str,
    output: str | Path,
    *,
    force_low_memory: bool = False,
    sort: bool = True,
) -> dict:
    """
    Convert a GTF file (local path or URL) to GFF3.

    Parameters
    ----------
    source
        Local path or http(s)/ftp URL (plain or gzip).
    output
        Output GFF3 path; use ``.gz`` suffix for gzip compression.
    force_low_memory
        When True, always sort by seqid into a temporary file before converting.
    sort
        When False (``--no-sort``), keep row order; skip seqid sort even for large files.

    Returns
    -------
    dict
        Summary with ``feature_count``, ``genes_with_cds``, ``transcripts_with_cds``, ``mode``.
    """
    source = str(source)
    output = Path(output)
    info = inspect_source(source)
    mode_label, _ = describe_convert_mode(
        source, force_low_memory=force_low_memory, sort=sort
    )

    needs_tmp = info.is_url or (sort and (force_low_memory or is_big(info)))

    if needs_tmp:
        with tempfile.TemporaryDirectory(prefix="gffy-convert-") as tmp_str:
            input_path = _prepare_input_path(
                source,
                info,
                Path(tmp_str),
                force_low_memory=force_low_memory,
                sort=sort,
            )
            feature_count, genes_with_cds, transcripts_with_cds = _run_conversion(
                input_path, output
            )
    else:
        input_path = Path(info.local_path or source)
        feature_count, genes_with_cds, transcripts_with_cds = _run_conversion(
            input_path, output
        )

    return {
        "feature_count": feature_count,
        "genes_with_cds": genes_with_cds,
        "transcripts_with_cds": transcripts_with_cds,
        "mode": mode_label,
    }
