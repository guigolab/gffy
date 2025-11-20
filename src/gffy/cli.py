#!/usr/bin/env python3
"""
Command-line interface for gffy - GFF3 Genomic Statistics Calculator
"""

import argparse
import json
import sys
from pathlib import Path

from gffy import compute_gff_stats


def main():
    parser = argparse.ArgumentParser(
        description="Compute comprehensive statistics from GFF3 genomic annotation files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s https://example.com/annotation.gff3.gz
  %(prog)s /path/to/local/annotation.gff3 --output stats.json
  %(prog)s annotation.gff3 --pretty --gzipped
        """
    )

    parser.add_argument(
        "gff_source",
        help="URL or local path to GFF3 file (may be compressed with .gz)"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output file path (default: stdout)"
    )

    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output with indentation"
    )

    parser.add_argument(
        "--gzipped",
        action="store_true",
        help="Gzip the input file"
    )

    args = parser.parse_args()

    try:
        # Redirect print statements from compute_gff_stats to stderr
        import builtins
        original_print = builtins.print

        def stderr_print(*args, **kwargs):
            kwargs.setdefault('file', sys.stderr)
            original_print(*args, **kwargs)

        builtins.print = stderr_print

        stats = compute_gff_stats(args.gff_source)

        # Restore original print
        builtins.print = original_print

        if not stats:
            print("Error: Failed to compute statistics", file=sys.stderr)
            sys.exit(1)

        if args.output:
            with open(args.output, 'w') as f:
                if args.pretty:
                    json.dump(stats, f, indent=2)
                else:
                    json.dump(stats, f)
            print(f"Statistics written to {args.output}")
        else:
            if args.pretty:
                print(json.dumps(stats, indent=2))
            else:
                print(json.dumps(stats))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
