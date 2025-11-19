from array import array
from collections import defaultdict
import gzip
import io
import os
from typing import Iterable, Optional

import requests

from .tools.helpers import (
    init_category_statistics_dictionary,
    init_orphan_feature,
    parse_gff_line_fast,
    process_feature,
)

SKIP_FEATURES = {"region", "chromosome", "scaffold"} #this speeds up the parsing process by skipping these features (can be around 500k lines saved)
CATEGORY_NAMES = ("coding", "long_non_coding", "short_non_coding", "pseudogene")


def _length_summary(
    values: array,
    decimals: int = 2,
) -> dict:
    """Return min/max/mean/median for an array using a single pass."""
    if not values:
        return {"min": 0, "max": 0, "mean": 0, "median": 0}

    iterator = iter(values)
    first_value = next(iterator)
    total = first_value
    min_value = first_value
    max_value = first_value
    collected = [first_value]

    for value in iterator:
        collected.append(value)
        total += value
        if value < min_value:
            min_value = value
        elif value > max_value:
            max_value = value

    collected.sort()
    count = len(collected)
    mid = count // 2
    if count % 2:
        median_value = collected[mid]
    else:
        median_value = (collected[mid - 1] + collected[mid]) / 2

    return {
        "min": int(min_value),
        "max": int(max_value),
        "mean": round(total / count, decimals),
        "median": round(median_value, decimals),
    }


def _init_category_store() -> dict:
    """Return a dictionary with initialized stats for every category."""
    return {name: init_category_statistics_dictionary() for name in CATEGORY_NAMES}


def _density(count: int, base: int) -> float:
    """Safe density helper that guards against division by zero."""
    return round(count / base, 2) if base else 0.0


def _build_feature_stats(
    count: int,
    denominator: int,
    lengths: array,
    concatenated_lengths: Optional[array] = None,
) -> dict:
    """Produce the standard feature stats structure used in the final JSON."""
    feature_stats = {
        "count": count,
        "density": _density(count, denominator),
        "length": _length_summary(lengths),
    }
    if concatenated_lengths is not None:
        feature_stats["length_concatenated"] = _length_summary(concatenated_lengths)
    return feature_stats


def _update_transcript_stats(transcript_info, transcript_stats: dict) -> None:
    """Append per-transcript measurements into the aggregate stats bucket."""
    transcript_stats["gene_ids"].add(transcript_info.gene_id)
    transcript_stats["length_list"].append(transcript_info.length)
    transcript_stats["exon_length_list"].extend(transcript_info.exons_lengths)
    transcript_stats["spliced_length_list"].append(transcript_info.exon_len_sum)

    if transcript_info.cds_count > 0:
        transcript_stats["cds_length_list"].extend(transcript_info.cds_lengths)
        transcript_stats["protein_length_list"].append(transcript_info.cds_len_sum // 3)

    if transcript_info.exon_count > 1:
        transcript_stats["intron_length_list"].append(
            transcript_info.length - transcript_info.exon_len_sum
        )


def compute_gff_stats(gff_source: str) -> dict:
    """
    Optimized GFF statistics computation from URL or local file path.

    Args:
        gff_source: URL (http/https/ftp) or local file path to GFF3 file (may be compressed)

    Returns:
        Dictionary containing comprehensive statistics for coding genes, long non-coding genes,
        short non-coding genes, and pseudogenes with transcript, exon, intron, and CDS metrics.

    Raises:
        FileNotFoundError: If local file path does not exist
        requests.RequestException: If URL cannot be accessed
        Exception: For other processing errors
    """
    print(f"Processing GFF from: {gff_source}")

    # Core accumulators populated during parsing
    roots = {}
    id_to_root = {}  # feature id -> root gene id
    transcripts = {}  # transcript id -> Transcript object

    # Deferred feature handling (child arrives before parent)
    waiting_by_parent: dict[str, list] = defaultdict(list)
    pending_orphans = []  # keep list for final reporting/debug output

    # Basic transport detection (explicit protocols avoid confusing local paths)
    is_remote_source = gff_source.startswith(("http://", "https://"))
    response = None
    file_obj = None

    if is_remote_source:
        print("Fetching from URL...")
        response = requests.get(gff_source, stream=True, timeout=120)
        response.raise_for_status()
        response.raw.decode_content = True
        file_obj = response.raw
    else:
        print("Reading local file...")
        if not os.path.isfile(gff_source):
            raise FileNotFoundError(f"GFF file not found: {gff_source}")
        file_obj = open(gff_source, "rb")

    # Auto-detect gzip compression (URLs rely on extension; local files inspect magic bytes)
    is_gzipped = False
    if is_remote_source:
        is_gzipped = gff_source.endswith(".gz") or gff_source.endswith(".gzip")
    else:
        magic_bytes = file_obj.read(2)
        is_gzipped = magic_bytes == b"\x1f\x8b"
        if hasattr(file_obj, "seek"):
            file_obj.seek(0)

    if is_gzipped:
        print("Detected gzip compression")
        gff_file_cm = gzip.GzipFile(fileobj=file_obj)
    else:
        print("Processing as uncompressed file")
        gff_file_cm = file_obj

    skip_features = SKIP_FEATURES
    parse_line = parse_gff_line_fast
    process = process_feature
    create_orphan = init_orphan_feature

    def handle_waiting_children(parent_id: Optional[str]) -> None:
        """Try to attach every deferred child whose parent just became known."""
        if not parent_id:
            return
        pending_children = waiting_by_parent.pop(parent_id, None)
        if not pending_children:
            return
        for orphan in pending_children:
            if orphan.resolved:
                continue
            processed_child = process(
                orphan.feature_id,
                orphan.feature_type,
                orphan.length,
                orphan.parent_ids,
                orphan.biotype,
                roots,
                id_to_root,
                transcripts,
            )
            if processed_child:
                orphan.resolved = True
                handle_waiting_children(orphan.feature_id)

    def iterate_stream(line_iterable: Iterable[str]) -> None:
        """Parse every row, attaching features or queuing them until their parent exists."""
        append_pending_orphan = pending_orphans.append
        waiting_lookup = waiting_by_parent
        handle_waiting = handle_waiting_children
        skip_feature_set = skip_features
        parse_columns = parse_line
        process_feature_local = process

        for line in line_iterable:
            if line.startswith("#"):
                continue

            # Only split the eight required columns; the ninth (attributes) stays intact
            cols = line.split("\t", 8)
            if len(cols) < 9:
                continue

            feature_type_raw = cols[2]
            if feature_type_raw in skip_feature_set:
                continue

            # Parse GFF attributes (ID, Parent, biotypes, etc.) only for relevant features
            feature_type, length, parent_ids, biotype, feature_id = parse_columns(cols)
            processed = process_feature_local(
                feature_id,
                feature_type,
                length,
                parent_ids,
                biotype,
                roots,
                id_to_root,
                transcripts,
            )

            if processed:
                # The feature is now anchored; unblock any descendants that were waiting
                handle_waiting(feature_id)
            elif parent_ids:
                # Parent hasn't been seen yet; store for deferred resolution
                orphan = create_orphan(feature_id, feature_type, length, parent_ids, biotype)
                append_pending_orphan(orphan)
                for parent_id in parent_ids:
                    waiting_lookup[parent_id].append(orphan)

    try:
        with gff_file_cm as stream:
            if isinstance(stream, io.TextIOBase):
                iterate_stream(stream)
            else:
                with io.TextIOWrapper(stream, encoding="utf-8", errors="ignore") as text_stream:
                    iterate_stream(text_stream)

        if pending_orphans:
            # Report any features that never found a parent, useful for debugging malformed GFFs
            unresolved = [orphan for orphan in pending_orphans if not orphan.resolved]
            resolved_count = len(pending_orphans) - len(unresolved)
            if resolved_count:
                print(f"Resolved {resolved_count} deferred features during parsing")
            if unresolved:
                print(f"Warning: {len(unresolved)} orphans could not be resolved")
                for orphan in unresolved[:5]:
                    print(
                        f"Orphan: {orphan.feature_id} {orphan.feature_type} "
                        f"{orphan.length} {orphan.parent_ids} {orphan.biotype}"
                    )


        # Set categories for all genes and lengths
        categories_dict = _init_category_store()

        for gene_info in roots.values():
            gene_info.set_category()
            if not gene_info.category:
                continue
            categories_dict[gene_info.category]["gene_length_list"].append(gene_info.length)
        
        # Process all transcripts to compute aggregate statistics
        get_root = roots.get
        for transcript_info in transcripts.values():
            root = get_root(transcript_info.gene_id)
            if not root or not root.category:
                continue

            transcript_stats = categories_dict[root.category]["types"][transcript_info.type]
            _update_transcript_stats(transcript_info, transcript_stats)
            
        # Build final results dictionary
        def build_category(category_name: str) -> dict:
            """Assemble the final metrics for a single gene category."""
            category_stats = categories_dict[category_name]
            gene_lengths = category_stats["gene_length_list"]

            category_obj = {
                "count": len(gene_lengths),
                "length": _length_summary(gene_lengths),
                "transcripts": {},
            }

            for transcript_type, transcript_stats in category_stats["types"].items():
                # Transcript stats step
                transcript_length_list = transcript_stats["length_list"]
                transcript_count = len(transcript_length_list)
                gene_id_count = len(transcript_stats["gene_ids"])
                transcript_density = _density(transcript_count, gene_id_count)
                type_entry = {
                    "count": transcript_count,
                    "density": transcript_density,
                    "length": _length_summary(transcript_length_list),
                    "features": {},
                }

                # Exon stats step
                exon_length_list = transcript_stats["exon_length_list"]
                exon_count = len(exon_length_list)
                type_entry["features"]["exon"] = _build_feature_stats(
                    exon_count,
                    transcript_count,
                    exon_length_list,
                    transcript_stats["spliced_length_list"],
                )

                # Intron stats step
                intron_length_list = transcript_stats["intron_length_list"]
                if intron_length_list:
                    intron_count = len(intron_length_list)
                    type_entry["features"]["intron"] = _build_feature_stats(
                        intron_count,
                        transcript_count,
                        intron_length_list,
                    )

                # CDS stats step
                cds_length_list = transcript_stats["cds_length_list"]
                if cds_length_list:
                    cds_count = len(cds_length_list)
                    type_entry["features"]["cds"] = _build_feature_stats(
                        cds_count,
                        transcript_count,
                        cds_length_list,
                        transcript_stats["protein_length_list"],
                    )

                category_obj["transcripts"][transcript_type] = type_entry

            return category_obj

        return {
            "coding_genes": build_category("coding"),
            "long_non_coding_genes": build_category("long_non_coding"),
            "short_non_coding_genes": build_category("short_non_coding"),
            "pseudogenes": build_category("pseudogene"),
        }

    except Exception as e:
        print(f"✗ Error: {str(e)}")
        import traceback

        traceback.print_exc()
        return {}
    finally:
        # Ensure network/file handles are released even if the parse fails
        if not is_remote_source and file_obj:
            try:
                file_obj.close()
            except Exception:
                pass
        if response is not None:
            response.close()
