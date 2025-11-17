import gzip
import os
import requests
from collections import defaultdict
import statistics
from .tools.helpers import (
    categorize_roots,
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
        Dictionary containing comprehensive statistics for coding genes, long non-coding genes,
        short non-coding genes, and pseudogenes with transcript, exon, intron, and CDS metrics.

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

        # Categorize gene roots into coding, long_non_coding, short_non_coding, and pseudogene categories
        root_to_category = categorize_roots(roots, transcripts)
        get_category = root_to_category.get

        # Precompute per-category gene lengths to avoid recomputation
        cat_gene_lengths = {
            "coding": [],
            "long_non_coding": [],
            "short_non_coding": [],
            "pseudogene": [],
        }
        for rid, info in roots.items():
            cat = get_category(rid)
            if cat:
                cat_gene_lengths[cat].append(info.length)

        # Initialize per-category statistics aggregates
        cat_type_stats = {
            "coding": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "transcript_lengths": [],
                "type_transcripts": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_exon_counts": defaultdict(list),
                "type_exon_lens": defaultdict(list),
                "type_concat_exon_lens": defaultdict(list),
                "type_intron_counts": defaultdict(list),
                "type_intron_lens": defaultdict(list),
                "type_concat_intron_lens": defaultdict(list),
                "type_cds_counts": defaultdict(list),
                "type_cds_lens": defaultdict(list),
                "type_concat_cds_lens": defaultdict(list),
            },
            "long_non_coding": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "transcript_lengths": [],
                "type_transcripts": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_exon_counts": defaultdict(list),
                "type_exon_lens": defaultdict(list),
                "type_concat_exon_lens": defaultdict(list),
                "type_intron_counts": defaultdict(list),
                "type_intron_lens": defaultdict(list),
                "type_concat_intron_lens": defaultdict(list),
            },
            "short_non_coding": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "transcript_lengths": [],
                "type_transcripts": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_exon_counts": defaultdict(list),
                "type_exon_lens": defaultdict(list),
                "type_concat_exon_lens": defaultdict(list),
                "type_intron_counts": defaultdict(list),
                "type_intron_lens": defaultdict(list),
                "type_concat_intron_lens": defaultdict(list),
            },
            "pseudogene": {
                "per_gene_type_counts": defaultdict(lambda: defaultdict(int)),
                "transcript_lengths": [],
                "type_transcripts": defaultdict(int),
                "type_span_lens": defaultdict(list),
                "type_exon_counts": defaultdict(list),
                "type_exon_lens": defaultdict(list),
                "type_concat_exon_lens": defaultdict(list),
                "type_intron_counts": defaultdict(list),
                "type_intron_lens": defaultdict(list),
                "type_concat_intron_lens": defaultdict(list),
            },
        }

        # Process all transcripts to compute aggregate statistics
        for tinfo in transcripts.values():
            category = get_category(tinfo.gene_id)
            if not category:
                continue

            cts = cat_type_stats[category]

            if tinfo.exons_flat and tinfo.type:
                ttype = tinfo.type
                exons_flat = tinfo.exons_flat
                exon_count = len(exons_flat) // 2
                
                cts["per_gene_type_counts"][tinfo.gene_id][ttype] += 1
                cts["type_transcripts"][ttype] += 1
                
                starts = exons_flat[0::2]
                ends = exons_flat[1::2]
                span_len = max(ends) - min(starts) + 1
                cts["type_span_lens"][ttype].append(span_len)
                cts["transcript_lengths"].append(span_len)
                
                # Exon statistics
                cts["type_exon_counts"][ttype].append(exon_count)
                cts["type_concat_exon_lens"][ttype].append(tinfo.exon_len_sum)
                
                for j in range(exon_count):
                    exon_len = ends[j] - starts[j] + 1
                    cts["type_exon_lens"][ttype].append(exon_len)
                
                # Intron statistics
                if exon_count > 1:
                    indices = list(range(0, len(exons_flat), 2))
                    indices.sort(key=lambda i: exons_flat[i])
                    
                    intron_count = 0
                    intron_len_sum = 0
                    intron_lens = []
                    for j in range(len(indices) - 1):
                        i = indices[j]
                        i_next = indices[j + 1]
                        intron_len = exons_flat[i_next] - exons_flat[i + 1] - 1
                        if intron_len > 0:
                            intron_count += 1
                            intron_len_sum += intron_len
                            intron_lens.append(intron_len)
                    
                    if intron_count > 0:
                        cts["type_intron_counts"][ttype].append(intron_count)
                        cts["type_concat_intron_lens"][ttype].append(intron_len_sum)
                        cts["type_intron_lens"][ttype].extend(intron_lens)
                
                # CDS statistics
                if tinfo.cds_segments > 0 and "type_cds_counts" in cts:
                    cts["type_cds_counts"][ttype].append(tinfo.cds_segments)
                    cts["type_concat_cds_lens"][ttype].append(tinfo.cds_total_len)
                    
                    for cds_len in tinfo.cds_lens:
                        cts["type_cds_lens"][ttype].append(cds_len)

        # Build final results dictionary
        def build_category(category_name: str) -> dict:
            gene_lengths = cat_gene_lengths[category_name]
            if not gene_lengths:
                return {}

            gene_count = len(gene_lengths)
            cts = cat_type_stats[category_name]
            per_gene_type_counts = cts["per_gene_type_counts"]
            type_transcripts = cts["type_transcripts"]
            type_span_lens = cts["type_span_lens"]
            transcript_lengths = cts["transcript_lengths"]
            
            type_exon_counts = cts["type_exon_counts"]
            type_exon_lens = cts["type_exon_lens"]
            type_concat_exon_lens = cts["type_concat_exon_lens"]
            type_intron_counts = cts["type_intron_counts"]
            type_intron_lens = cts["type_intron_lens"]
            type_concat_intron_lens = cts["type_concat_intron_lens"]
            type_cds_counts = cts.get("type_cds_counts", defaultdict(list))
            type_cds_lens = cts.get("type_cds_lens", defaultdict(list))
            type_concat_cds_lens = cts.get("type_concat_cds_lens", defaultdict(list))

            child_totals = defaultdict(int)
            child_counts_per_gene = defaultdict(list)
            for type_counts in per_gene_type_counts.values():
                for ctype, cnt in type_counts.items():
                    child_totals[ctype] += cnt
                    child_counts_per_gene[ctype].append(cnt)

            total_transcripts = sum(type_transcripts.values())

            category_obj = {
                "count": gene_count,
                "length": {
                    "min": min(gene_lengths) if gene_lengths else 0,
                    "max": max(gene_lengths) if gene_lengths else 0,
                    "mean": round(sum(gene_lengths) / len(gene_lengths), 2) if gene_lengths else 0,
                    "median": statistics.median(gene_lengths) if gene_lengths else 0,
                },
                "transcripts": {
                    "count": total_transcripts,
                    "density": round(total_transcripts / gene_count, 2) if gene_count else 0,
                    "length": {
                        "min": min(transcript_lengths) if transcript_lengths else 0,
                        "max": max(transcript_lengths) if transcript_lengths else 0,
                        "mean": round(statistics.mean(transcript_lengths), 2) if transcript_lengths else 0,
                        "median": statistics.median(transcript_lengths) if transcript_lengths else 0,
                    },
                    "by_type": {},
                },
            }
            
            for ctype, total in child_totals.items():
                counts = child_counts_per_gene[ctype]
                type_trans = type_transcripts.get(ctype, 0)

                span_lens_list = type_span_lens.get(ctype, [])
                
                type_entry = {
                    "count": total,
                    "density": round(sum(counts) / len(counts), 2) if counts else 0,
                    "length": {
                        "min": min(span_lens_list) if span_lens_list else 0,
                        "max": max(span_lens_list) if span_lens_list else 0,
                        "mean": round(statistics.mean(span_lens_list), 2) if span_lens_list else 0,
                        "median": statistics.median(span_lens_list) if span_lens_list else 0,
                    },
                    "features": {},
                }

                # Exon statistics
                exon_counts_list = type_exon_counts.get(ctype, [])
                exon_lens_list = type_exon_lens.get(ctype, [])
                concat_exon_lens_list = type_concat_exon_lens.get(ctype, [])
                
                if exon_counts_list:
                    total_exons = sum(exon_counts_list)
                    type_entry["features"]["exon"] = {
                        "count": total_exons,
                        "density": round(statistics.mean(exon_counts_list), 2),
                        "length": {
                            "min": min(exon_lens_list) if exon_lens_list else 0,
                            "max": max(exon_lens_list) if exon_lens_list else 0,
                            "mean": round(statistics.mean(exon_lens_list), 2) if exon_lens_list else 0,
                            "median": round(statistics.median(exon_lens_list), 2) if exon_lens_list else 0,
                        },
                        "length_concatenated": {
                            "min": min(concat_exon_lens_list) if concat_exon_lens_list else 0,
                            "max": max(concat_exon_lens_list) if concat_exon_lens_list else 0,
                            "mean": round(statistics.mean(concat_exon_lens_list), 2) if concat_exon_lens_list else 0,
                            "median": statistics.median(concat_exon_lens_list) if concat_exon_lens_list else 0,
                        }
                    }
                
                # Intron statistics
                intron_counts_list = type_intron_counts.get(ctype, [])
                intron_lens_list = type_intron_lens.get(ctype, [])
                concat_intron_lens_list = type_concat_intron_lens.get(ctype, [])
                
                if intron_counts_list:
                    total_introns = sum(intron_counts_list)
                    type_entry["features"]["intron"] = {
                        "count": total_introns,
                        "density": round(statistics.mean(intron_counts_list), 2),
                        "length": {
                            "min": min(intron_lens_list) if intron_lens_list else 0,
                            "max": max(intron_lens_list) if intron_lens_list else 0,
                            "mean": round(statistics.mean(intron_lens_list), 2) if intron_lens_list else 0,
                            "median": round(statistics.median(intron_lens_list), 2) if intron_lens_list else 0,
                        },
                        "length_concatenated": {
                            "min": min(concat_intron_lens_list) if concat_intron_lens_list else 0,
                            "max": max(concat_intron_lens_list) if concat_intron_lens_list else 0,
                            "mean": round(statistics.mean(concat_intron_lens_list), 2) if concat_intron_lens_list else 0,
                            "median": statistics.median(concat_intron_lens_list) if concat_intron_lens_list else 0,
                        }
                    }
                
                # CDS statistics
                cds_counts_list = type_cds_counts.get(ctype, [])
                cds_lens_list = type_cds_lens.get(ctype, [])
                concat_cds_lens_list = type_concat_cds_lens.get(ctype, [])
                
                if cds_counts_list:
                    total_cds = sum(cds_counts_list)
                    type_entry["features"]["cds"] = {
                        "count": total_cds,
                        "density": round(statistics.mean(cds_counts_list), 2),
                        "length": {
                            "min": min(cds_lens_list) if cds_lens_list else 0,
                            "max": max(cds_lens_list) if cds_lens_list else 0,
                            "mean": round(statistics.mean(cds_lens_list), 2) if cds_lens_list else 0,
                            "median": round(statistics.median(cds_lens_list), 2) if cds_lens_list else 0,
                        },
                        "length_concatenated": {
                            "min": min(concat_cds_lens_list) if concat_cds_lens_list else 0,
                            "max": max(concat_cds_lens_list) if concat_cds_lens_list else 0,
                            "mean": round(statistics.mean(concat_cds_lens_list), 2) if concat_cds_lens_list else 0,
                            "median": statistics.median(concat_cds_lens_list) if concat_cds_lens_list else 0,
                        }
                    }

                category_obj["transcripts"]["by_type"][ctype] = type_entry

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
