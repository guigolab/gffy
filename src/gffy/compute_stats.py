from array import array
import statistics
from .tools.helpers import (
    init_category_statistics_dictionary,
    init_orphan_feature,
    parse_gff_line_fast,
    process_feature,
    resolve_orphans,
)
import requests
import os
import gzip

SKIP_FEATURES = {"region", "chromosome", "scaffold"} #this speeds up the parsing process by skipping these features (can be around 500k lines saved)


def _length_summary(
    values: array,
    decimals: int = 2,
) -> dict:
    if not values:
        return {"min": 0, "max": 0, "mean": 0, "median": 0}

    mean = round(statistics.mean(values), decimals)
    median_value = round(statistics.median(values), decimals)

    return {
        "min": int(min(values)),
        "max": int(max(values)),
        "mean": mean,
        "median": median_value,
    }


def compute_gff_stats(gff_source: str) -> dict:
    """
    Optimized GFF statistics computation from URL or local file path.

    Args:
        gff_source: URL (http/https/ftp) or local file path to GFF3 file (may be compressed)
        skip_features: optional set of feature types to ignore during parsing
        is_gzipped: when True, skip gzip autodetection (mainly for callers that already know)

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
    id_to_root = {} #feature id to root id mapping
    transcripts = {} #transcript dictionary -> key is transcript id, value is @Transcript object
    orphans = [] #list of orphan features

    #resolve gff path
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

    try:
        # Process GFF file line by line using context manager when possible
        with gff_file_cm as stream:
            for line in stream:
                line = line.decode("utf-8", errors="ignore")
                if line.startswith("#"):
                    continue

                cols = line.split("\t")
                if len(cols) < 9:
                    continue
                
                # Skip region features (using set for O(1) lookup)
                if cols[2] in SKIP_FEATURES:
                    continue

                feature_type, length, parent_ids, biotype, feature_id = parse_gff_line_fast(
                    cols
                )

                processed = process_feature(
                    feature_id,
                    feature_type,
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
                            feature_id, feature_type, length, parent_ids, biotype
                        )
                    )

        if file_needs_closing:
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
                #print first 5 orphans
                for orphan in orphans[:5]:
                    print(f"Orphan: {orphan.feature_id} {orphan.feature_type} {orphan.length} {orphan.parent_ids} {orphan.biotype}")


        # Set categories for all genes and lenghts
        categories_dict = {
            "coding": init_category_statistics_dictionary(),
            "long_non_coding": init_category_statistics_dictionary(),
            "short_non_coding": init_category_statistics_dictionary(),
            "pseudogene": init_category_statistics_dictionary(),
        }

        for gene_info in roots.values():
            gene_info.set_category()
            if not gene_info.category:
                continue
            categories_dict[gene_info.category]["gene_length_list"].append(gene_info.length)
        
        # Process all transcripts to compute aggregate statistics
        get_root = roots.get
        for transcript_info in transcripts.values():
            category = get_root(transcript_info.gene_id).category
            if not category:
                continue

            transcript_types_stats = categories_dict[category]['types'][transcript_info.type]
            transcript_types_stats["gene_ids"].add(transcript_info.gene_id)
            transcript_types_stats["length_list"].append(transcript_info.length)
            transcript_types_stats["exon_length_list"].extend(transcript_info.exons_lengths)
            transcript_types_stats["spliced_length_list"].append(transcript_info.exon_len_sum)
            
            if transcript_info.cds_count > 0:
                transcript_types_stats["cds_length_list"].extend(transcript_info.cds_lengths)
                transcript_types_stats["protein_length_list"].append((transcript_info.cds_len_sum // 3))
            if transcript_info.exon_count > 1:
                transcript_types_stats["intron_length_list"].append(transcript_info.length - transcript_info.exon_len_sum)
            
        # Build final results dictionary
        def build_category(category_name: str) -> dict:
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
                type_entry = {
                    "count": len(transcript_length_list),
                    "density": round(len(transcript_length_list) / len(transcript_stats['gene_ids']), 2),
                    "length": _length_summary(transcript_length_list),
                    "features": {},
                }

                # Exon stats step
                exon_length_list = transcript_stats["exon_length_list"]
                type_entry["features"]["exon"] = {
                    "count": len(exon_length_list),
                    "density": round(len(exon_length_list) / type_entry["count"], 2),
                    "length": _length_summary(exon_length_list),
                    "length_concatenated": _length_summary(transcript_stats["spliced_length_list"]),
                }

                # Intron stats step
                intron_length_list = transcript_stats["intron_length_list"]
                if intron_length_list:
                    type_entry["features"]["intron"] = {
                        "count": len(intron_length_list),
                        "density": round(len(intron_length_list) / type_entry["count"], 2),
                        "length": _length_summary(intron_length_list),
                    }

                # CDS stats step
                cds_length_list = transcript_stats["cds_length_list"]
                if cds_length_list:
                    type_entry["features"]["cds"] = {
                        "count": len(cds_length_list),
                        "density": round(len(cds_length_list) / type_entry["count"], 2),
                        "length": _length_summary(cds_length_list),
                        "length_concatenated": _length_summary(transcript_stats["protein_length_list"]),
                    }

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
