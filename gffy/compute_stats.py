import gzip
import os
import requests
from collections import defaultdict
import statistics
from .tools.helpers import (
    categorize_roots,
    init_category_totals,
    init_orphan_feature,
    parse_gff_line_fast,
    process_feature,
    resolve_orphans,
)

SKIP_FEATURES = {"region", "chromosome", "scaffold"}


def compute_gff_stats(gff_source: str) -> dict:
    """
    Optimized GFF statistics computation from URL or local file path.

    Args:
        gff_source: URL (http/https/ftp) or local file path to GFF3 file (may be compressed)

    Returns:
        Dictionary containing comprehensive statistics for coding genes, non-coding genes,
        and pseudogenes with transcript, exon, intron, and CDS metrics.

    Raises:
        FileNotFoundError: If local file path does not exist
        requests.RequestException: If URL cannot be accessed
        Exception: For other processing errors
    """
    print(f"Processing GFF from: {gff_source}")

    roots = {}
    id_to_root = {}
    transcripts = {}
    orphans = []

    # Pre-computed skip set for feature types to ignore
    skip_features = SKIP_FEATURES

    try:
        # Determine if source is URL or local file and open appropriately
        if gff_source.startswith(("http://", "https://", "ftp://")):
            # Handle remote URL
            print("Fetching from URL...")
            response = requests.get(gff_source, stream=True, timeout=120)
            response.raise_for_status()
            response.raw.decode_content = True
            file_obj = response.raw
            is_remote = True
        else:
            # Handle local file
            print("Reading local file...")
            if not os.path.isfile(gff_source):
                raise FileNotFoundError(f"GFF file not found: {gff_source}")
            file_obj = open(gff_source, "rb")
            is_remote = False

        # Auto-detect gzip compression
        is_gzipped = False
        if is_remote:
            # For remote files, check if URL ends with .gz or assume gzipped for common extensions
            is_gzipped = gff_source.endswith(".gz") or gff_source.endswith(".gzip")
        else:
            # For local files, check magic bytes
            if hasattr(file_obj, "read"):
                # Peek at first two bytes to check for gzip magic number
                magic_bytes = file_obj.read(2)
                is_gzipped = magic_bytes == b"\x1f\x8b"
                # Reset file pointer
                if hasattr(file_obj, "seek"):
                    file_obj.seek(0)

        # Open file with appropriate compression handling
        if is_gzipped:
            print("Detected gzip compression")
            gff_file_cm = gzip.GzipFile(fileobj=file_obj)
        else:
            print("Processing as uncompressed file")
            gff_file_cm = file_obj

        # Track if we need to manually close the file object (for local uncompressed files)
        file_needs_closing = not is_remote and not is_gzipped

        # Process GFF file line by line using context manager when possible
        with gff_file_cm as gff_file:
            for line_bytes in gff_file:
                # Quick comment check
                line = line_bytes.decode("utf-8", errors="ignore")
                if not line or line[0] == "#":
                    continue

                line = line.strip()
                if not line:
                    continue

                cols = line.split("\t")
                if len(cols) < 9:
                    continue

                # Skip region features (using set for O(1) lookup)
                if cols[2] in skip_features:
                    continue

                feature_type, start, end, parent_ids, biotype, feature_id = parse_gff_line_fast(
                    cols
                )
                length = end - start + 1

                processed = process_feature(
                    feature_id,
                    feature_type,
                    start,
                    end,
                    length,
                    parent_ids,
                    biotype,
                    roots,
                    id_to_root,
                    transcripts,
                )

                if not processed and parent_ids:
                    orphans.append(
                        init_orphan_feature(
                            feature_id, feature_type, start, end, parent_ids, biotype
                        )
                    )

        # Clean up file resources for local files that weren't handled by context manager
        if file_needs_closing and hasattr(file_obj, "close"):
            file_obj.close()

        # Resolve orphan features (features that appeared before their parents)
        if orphans:
            print(f"Initial orphans: {len(orphans)}")
            max_iterations = 20
            iteration = 0

            while orphans and iteration < max_iterations:
                iteration += 1
                prev_count = len(orphans)
                orphans = resolve_orphans(orphans, roots, id_to_root, transcripts)
                resolved = prev_count - len(orphans)

                if resolved > 0:
                    print(f"Iter {iteration}: {resolved} orphans resolved, {len(orphans)} left")

                if prev_count == len(orphans):
                    break

            if orphans:
                print(f"Warning: {len(orphans)} orphans could not be resolved")

        # Categorize gene roots into coding, non-coding, and pseudogene categories
        root_to_category = categorize_roots(roots)
        get_category = root_to_category.get

        # Initialize totals for each gene category
        category_totals = init_category_totals()

        # Precompute per-category gene lengths to avoid recomputation
        cat_gene_lengths = {
            "coding": [],
            "non_coding": [],
            "pseudogene": [],
        }
        for rid, info in roots.items():
            cat = get_category(rid)
            if cat:
                cat_gene_lengths[cat].append(info.length)

        # Initialize per-category per-type statistics aggregates
        # This avoids rescanning transcripts during final statistics building
        cat_type_stats = {
            "coding": {
                "per_gene_type_counts": defaultdict(
                    lambda: defaultdict(int)
                ),  # gene_id -> type -> count
                "type_transcripts": defaultdict(int),  # type -> total transcripts
                "type_exons": defaultdict(int),  # type -> total exons
                "type_span_sum": defaultdict(int),  # type -> sum of genomic spans
                "type_spliced_sum": defaultdict(int),  # type -> sum of spliced lengths
                "type_exon_len_sum": defaultdict(int),  # type -> sum of exon lengths
                "type_span_lens": defaultdict(list),  # type -> list of span lengths
                "type_spliced_lens": defaultdict(list),  # type -> list of spliced lengths
                "type_exon_lens": defaultdict(list),  # type -> list of exon lengths
            },
            "non_coding": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "type_transcripts": defaultdict(int),
                "type_exons": defaultdict(int),
                "type_span_sum": defaultdict(int),
                "type_spliced_sum": defaultdict(int),
                "type_exon_len_sum": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_spliced_lens": defaultdict(list),
                "type_exon_lens": defaultdict(list),
            },
            "pseudogene": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "type_transcripts": defaultdict(int),
                "type_exons": defaultdict(int),
                "type_span_sum": defaultdict(int),
                "type_spliced_sum": defaultdict(int),
                "type_exon_len_sum": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_spliced_lens": defaultdict(list),
                "type_exon_lens": defaultdict(list),
            },
        }

        # Process all transcripts to compute aggregate statistics
        for tinfo in transcripts.values():
            category = get_category(tinfo.gene_id)
            if not category:
                continue

            cat_totals = category_totals[category]
            cts = cat_type_stats[category]

            if tinfo.cds_total_len > 0:
                cat_totals["cds"] += tinfo.cds_segments
                cat_totals["cds_len_sum"] += tinfo.cds_total_len
                cat_totals["cds_transcripts"] += 1

            exons_flat = tinfo.exons_flat
            exon_count = len(exons_flat) // 2
            cat_totals["exons"] += exon_count
            cat_totals["exon_len_sum"] += tinfo.exon_len_sum

            # Optimized intron calculation
            if exon_count > 1:
                indices = list(range(0, len(exons_flat), 2))
                indices.sort(key=lambda i: exons_flat[i])

                for j in range(len(indices) - 1):
                    i = indices[j]
                    i_next = indices[j + 1]
                    intron_len = exons_flat[i_next] - exons_flat[i + 1] - 1
                    if intron_len > 0:
                        cat_totals["introns"] += 1
                        cat_totals["intron_len_sum"] += intron_len
                        # Store intron length for median later
                        intron_lens = cat_totals.setdefault("intron_lens", [])
                        intron_lens.append(intron_len)

            # Store exon lengths for median
            if exon_count:
                exon_lens = cat_totals.setdefault("exon_lens", [])
                # push each exon length quickly
                for k in range(0, len(exons_flat), 2):
                    exon_lens.append(exons_flat[k + 1] - exons_flat[k] + 1)

            # Store CDS lengths for median (per segment)
            if tinfo.cds_segments:
                cds_lens = cat_totals.setdefault("cds_lens", [])
                cds_lens.extend(tinfo.cds_lens)

            # Per-type aggregates (only for transcripts with exons and with a defined type)
            if tinfo.exons_flat and tinfo.type:
                ttype = tinfo.type
                exons_flat = tinfo.exons_flat
                exon_count = len(exons_flat) // 2
                # counts per gene for per_gene ratios
                cts["per_gene_type_counts"][tinfo.gene_id][ttype] += 1
                # totals
                cts["type_transcripts"][ttype] += 1
                cts["type_exons"][ttype] += exon_count
                # span and spliced
                starts = exons_flat[0::2]
                ends = exons_flat[1::2]
                span_len = max(ends) - min(starts) + 1
                cts["type_span_sum"][ttype] += span_len
                cts["type_spliced_sum"][ttype] += tinfo.exon_len_sum
                cts["type_span_lens"][ttype].append(span_len)
                cts["type_spliced_lens"][ttype].append(tinfo.exon_len_sum)
                # exon lengths
                for k in range(0, len(exons_flat), 2):
                    elen = exons_flat[k + 1] - exons_flat[k] + 1
                    cts["type_exon_len_sum"][ttype] += elen
                    cts["type_exon_lens"][ttype].append(elen)

        # Build final results dictionary
        def build_category(category_name: str) -> dict:
            gene_lengths = cat_gene_lengths[category_name]
            if not gene_lengths:
                return {}

            gene_count = len(gene_lengths)
            cts = cat_type_stats[category_name]
            per_gene_type_counts = cts["per_gene_type_counts"]
            type_transcripts = cts["type_transcripts"]
            type_exons = cts["type_exons"]
            type_span_sum = cts["type_span_sum"]
            type_spliced_sum = cts["type_spliced_sum"]
            type_exon_len_sum = cts["type_exon_len_sum"]
            type_span_lens = cts["type_span_lens"]
            type_spliced_lens = cts["type_spliced_lens"]
            type_exon_lens = cts["type_exon_lens"]

            child_totals = defaultdict(int)
            child_counts_per_gene = defaultdict(list)
            for type_counts in per_gene_type_counts.values():
                for ctype, cnt in type_counts.items():
                    child_totals[ctype] += cnt
                    child_counts_per_gene[ctype].append(cnt)

            total_transcripts = sum(type_transcripts.values())
            cat_totals = category_totals[category_name]
            get_lens = cat_totals.get
            # Medians
            gene_median = statistics.median(gene_lengths) if gene_lengths else 0
            exon_median = (
                statistics.median(get_lens("exon_lens", [])) if get_lens("exon_lens") else 0
            )
            intron_median = (
                statistics.median(get_lens("intron_lens", [])) if get_lens("intron_lens") else 0
            )
            cds_median = statistics.median(get_lens("cds_lens", [])) if get_lens("cds_lens") else 0

            category_obj = {
                "count": gene_count,
                "length_stats": {
                    "min": min(gene_lengths) if gene_lengths else 0,
                    "max": max(gene_lengths) if gene_lengths else 0,
                    "mean": round(sum(gene_lengths) / len(gene_lengths), 2) if gene_lengths else 0,
                    "median": gene_median,
                },
                "transcripts": {
                    "count": total_transcripts,
                    "per_gene": round(total_transcripts / gene_count, 2) if gene_count else 0,
                    "types": {},
                },
                "features": {
                    "exons": {
                        "count": cat_totals["exons"],
                        "length_stats": {
                            "mean": (
                                round(cat_totals["exon_len_sum"] / cat_totals["exons"], 2)
                                if cat_totals["exons"]
                                else 0
                            ),
                            "median": exon_median,
                        },
                    },
                    "introns": {
                        "count": cat_totals["introns"],
                        "length_stats": {
                            "mean": (
                                round(cat_totals["intron_len_sum"] / cat_totals["introns"], 2)
                                if cat_totals["introns"]
                                else 0
                            ),
                            "median": intron_median,
                        },
                    },
                },
            }
            if cat_totals["cds"] > 0:
                category_obj["features"]["cds"] = {
                    "count": cat_totals["cds"],
                    "length_stats": {
                        "mean": (
                            round(cat_totals["cds_len_sum"] / cat_totals["cds_transcripts"], 2)
                            if cat_totals["cds_transcripts"]
                            else 0
                        ),
                        "median": cds_median,
                    },
                }
            get_type_trans = type_transcripts.get
            for ctype, total in child_totals.items():
                counts = child_counts_per_gene[ctype]
                type_trans = get_type_trans(ctype, 0)

                exons_total_for_type = type_exons.get(ctype, 0)
                one_exon_per_transcript = type_trans > 0 and exons_total_for_type == type_trans

                type_entry = {
                    "count": total,
                    "per_gene": round(sum(counts) / len(counts), 2) if counts else 0,
                    "exons_per_transcript": (
                        round(exons_total_for_type / type_trans, 2) if type_trans else 0
                    ),
                    "length_stats": {
                        "mean": round(type_span_sum[ctype] / type_trans, 2) if type_trans else 0,
                        "median": (
                            statistics.median(type_span_lens[ctype])
                            if type_span_lens.get(ctype)
                            else 0
                        ),
                    },
                }

                if not one_exon_per_transcript:
                    type_entry["spliced_length_stats"] = {
                        "mean": round(type_spliced_sum[ctype] / type_trans, 2) if type_trans else 0,
                        "median": (
                            statistics.median(type_spliced_lens[ctype])
                            if type_spliced_lens.get(ctype)
                            else 0
                        ),
                    }
                    type_entry["exon_length_stats"] = {
                        "mean": (
                            round(type_exon_len_sum[ctype] / exons_total_for_type, 2)
                            if exons_total_for_type
                            else 0
                        ),
                        "median": (
                            statistics.median(type_exon_lens[ctype])
                            if type_exon_lens.get(ctype)
                            else 0
                        ),
                    }

                category_obj["transcripts"]["types"][ctype] = type_entry

            return category_obj

        return {
            "coding_genes": build_category("coding"),
            "non_coding_genes": build_category("non_coding"),
            "pseudogenes": build_category("pseudogene"),
        }

    except Exception as e:
        print(f"✗ Error: {str(e)}")
        import traceback

        traceback.print_exc()
        return {}
