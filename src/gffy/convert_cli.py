#!/usr/bin/env python3
"""Command-line interface for gffy-convert - GTF to GFF3 converter."""

import argparse
import sys
from pathlib import Path

from gffy.convert import convert_gtf_to_gff3, describe_convert_mode
from gffy._cli_common import add_usage_argument, print_resource_usage
from gffy.errors import print_classified_error


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="gffy-convert",
        description="Convert GTF annotations to GFF3 (ID/Parent attributes, biotype inference)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s annotation.gtf -o annotation.gff3
  %(prog)s annotation.gtf.gz -o annotation.gff3.gz
  %(prog)s https://example.com/annotation.gtf.gz -o local.gff3 --low-memory
  %(prog)s large.gtf.gz -o out.gff3 --low-memory
  %(prog)s annotation.gtf -o out.gff3 --no-sort
        """,
    )

    parser.add_argument(
        "gtf_source",
        help="URL or local path to GTF file (may be compressed with .gz)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output GFF3 file path (use .gz for gzip compression)",
    )
    parser.add_argument(
        "--low-memory",
        action="store_true",
        help="Force sort-to-temp by seqid before converting (large or unordered files)",
    )
    parser.add_argument(
        "--no-sort",
        action="store_true",
        help="Keep input row order (skip seqid sort; still stages URLs once)",
    )
    add_usage_argument(parser)

    args = parser.parse_args()

    try:
        mode_label, _info = describe_convert_mode(
            args.gtf_source,
            force_low_memory=args.low_memory,
            sort=not args.no_sort,
        )
        print(f"[gffy-convert] Mode: {mode_label}", file=sys.stderr)

        summary = convert_gtf_to_gff3(
            args.gtf_source,
            args.output,
            force_low_memory=args.low_memory,
            sort=not args.no_sort,
        )

        print(
            f"[gffy-convert] Wrote {summary['feature_count']} features to {args.output}",
            file=sys.stderr,
        )
        print(
            f"[gffy-convert] CDS-bearing: "
            f"{summary['transcripts_with_cds']} transcripts, "
            f"{summary['genes_with_cds']} genes",
            file=sys.stderr,
        )

    except FileNotFoundError as e:
        print_classified_error("[gffy-convert]", e)
        sys.exit(1)
    except Exception as e:
        print_classified_error("[gffy-convert]", e)
        sys.exit(1)
    finally:
        if args.usage:
            print_resource_usage("[gffy-convert]")


if __name__ == "__main__":
    main()
