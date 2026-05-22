#!/usr/bin/env python3
"""Command-line interface for gffy - GFF3 Genomic Statistics Calculator."""

import argparse
import json
import os
import sys
from pathlib import Path

from gffy import compute_gff_stats
from gffy.annotrieve_export import build_custom_annotation, write_annotrieve_json
from gffy.bulk import (
    compute_bulk_annotrieve_export,
    compute_bulk_stats,
    read_source_list,
    warn_source_list_issues,
)
from gffy._cli_common import add_usage_argument, print_resource_usage
from gffy.errors import print_classified_error
from gffy.stats import describe_compute_mode


def _bulk_progress(source: str, record: dict) -> None:
    elapsed = record.get("elapsed_ms", 0)
    if "error" in record:
        category = record.get("error_category", "unknown")
        detail = record.get("error", "")
        print(
            f"[gffy] err {category} {source} ({elapsed} ms): {detail}",
            file=sys.stderr,
        )
    else:
        print(f"[gffy] ok {source} ({elapsed} ms)", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="gffy",
        description="Compute GFF3 feature statistics (annotrieve GFFStats schema)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s https://example.com/annotation.gff3.gz
  %(prog)s /path/to/local/annotation.gff3 --output stats.json
  %(prog)s annotation.gff3 --pretty
  %(prog)s large.gff.gz --low-memory
  %(prog)s --from-file urls.txt -o stats.jsonl --workers 4
  %(prog)s ann.gff3 --annotrieve-json favorites.json
  %(prog)s --from-file urls.txt --annotrieve-jsonl import.jsonl
  cat urls.txt | %(prog)s --from-file - -o stats.jsonl
        """,
    )

    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "gff_source",
        nargs="?",
        help="URL or local path to GFF3 file (may be compressed with .gz)",
    )
    source_group.add_argument(
        "--from-file",
        metavar="PATH",
        help="Newline-separated list of paths/URLs (- for stdin); writes JSONL",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output file path (required for --from-file; default stdout otherwise)",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output with indentation (single-source only)",
    )
    parser.add_argument(
        "--low-memory",
        action="store_true",
        help="Force sort-to-temp + per-seqid mode (for large or unordered files)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Parallel workers for --from-file (default: 1; "
            f"max {min(os.cpu_count() or 1, 8)}). Each worker is a separate process."
        ),
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop bulk run on first failing source (default: record error and continue)",
    )
    parser.add_argument(
        "--failed-urls-out",
        type=Path,
        metavar="PATH",
        help=(
            "Write unfetchable URLs from bulk run to PATH "
            "(default: <output>.failed_urls.txt)"
        ),
    )
    parser.add_argument(
        "--annotrieve-json",
        type=Path,
        metavar="PATH",
        help="Write Annotrieve Favorites import JSON (single source)",
    )
    parser.add_argument(
        "--custom-name",
        metavar="NAME",
        help="Override custom_name for single-source --annotrieve-json",
    )
    parser.add_argument(
        "--annotrieve-jsonl",
        type=Path,
        metavar="PATH",
        help="Write import-ready JSONL for --from-file bulk (one record per line)",
    )
    add_usage_argument(parser)

    args = parser.parse_args()

    if args.custom_name and not args.annotrieve_json:
        print(
            "Error: --custom-name requires --annotrieve-json",
            file=sys.stderr,
        )
        sys.exit(2)
    if args.annotrieve_jsonl and not args.from_file:
        print(
            "Error: --annotrieve-jsonl requires --from-file",
            file=sys.stderr,
        )
        sys.exit(2)

    try:
        if args.from_file:
            if not args.output and not args.annotrieve_jsonl:
                print(
                    "Error: --output or --annotrieve-jsonl is required with --from-file",
                    file=sys.stderr,
                )
                sys.exit(2)
            if args.pretty:
                print(
                    "Error: --pretty is not supported with --from-file (JSONL output)",
                    file=sys.stderr,
                )
                sys.exit(2)

            list_result = read_source_list(args.from_file)
            warn_source_list_issues(list_result)
            print(
                f"[gffy] Bulk: {len(list_result.sources)} source(s), "
                f"workers={args.workers}",
                file=sys.stderr,
            )

            if args.output:
                summary = compute_bulk_stats(
                    list_result.sources,
                    args.output,
                    workers=args.workers,
                    force_low_memory=args.low_memory,
                    continue_on_error=not args.fail_fast,
                    failed_urls_output=args.failed_urls_out,
                    progress=_bulk_progress,
                    duplicates_skipped=len(list_result.duplicates),
                )
                print(
                    f"[gffy] Wrote {summary['total']} records to {args.output} "
                    f"({summary['succeeded']} ok, {summary['failed']} failed)",
                    file=sys.stderr,
                )
                if summary.get("duplicates_skipped"):
                    print(
                        f"[gffy] Skipped {summary['duplicates_skipped']} duplicate line(s)",
                        file=sys.stderr,
                    )
                if summary["failed"] > 0 and args.fail_fast:
                    sys.exit(1)

            if args.annotrieve_jsonl:
                ann_summary = compute_bulk_annotrieve_export(
                    list_result.sources,
                    args.annotrieve_jsonl,
                    workers=args.workers,
                    force_low_memory=args.low_memory,
                    continue_on_error=not args.fail_fast,
                    progress=_bulk_progress,
                )
                print(
                    f"[gffy] Wrote {ann_summary['succeeded']} Annotrieve record(s) to "
                    f"{args.annotrieve_jsonl} "
                    f"({ann_summary['failed']} failed, {ann_summary['total']} total)",
                    file=sys.stderr,
                )
                if ann_summary["failed"] > 0 and args.fail_fast:
                    sys.exit(1)
            return

        want_stats = bool(args.output) or not args.annotrieve_json

        if want_stats:
            mode_label, _info = describe_compute_mode(
                args.gff_source, force_low_memory=args.low_memory
            )
            print(f"[gffy] Mode: {mode_label}", file=sys.stderr)

            stats = compute_gff_stats(
                args.gff_source,
                force_low_memory=args.low_memory,
            )

            if args.output:
                with open(args.output, "w", encoding="utf-8") as f:
                    if args.pretty:
                        json.dump(stats, f, indent=2)
                    else:
                        json.dump(stats, f)
                print(f"Statistics written to {args.output}", file=sys.stderr)
            else:
                if args.pretty:
                    print(json.dumps(stats, indent=2))
                else:
                    print(json.dumps(stats))

        if args.annotrieve_json:
            record = build_custom_annotation(
                args.gff_source,
                custom_name_override=args.custom_name,
                force_low_memory=args.low_memory,
            )
            write_annotrieve_json(record, args.annotrieve_json)
            print(
                f"Annotrieve import JSON written to {args.annotrieve_json}",
                file=sys.stderr,
            )

    except FileNotFoundError as e:
        print_classified_error("[gffy]", e)
        sys.exit(1)
    except Exception as e:
        print_classified_error("[gffy]", e)
        sys.exit(1)
    finally:
        if args.usage:
            print_resource_usage("[gffy]")


if __name__ == "__main__":
    main()
