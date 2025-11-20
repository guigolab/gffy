#!/usr/bin/env python3
"""
Utility to compute statistics for a fixed GFF source and validate against schema.json.

This is primarily used by CI but can also be executed locally:

    python scripts/compute_and_validate.py
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

from jsonschema import Draft7Validator

from gffy import compute_gff_stats

GFF_SOURCE = (
    "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/"
    "GCA_000001405.29/ensembl/geneset/2024_11/genes.gff3.gz"
)
SCHEMA_PATH = Path("schema.json")
OUTPUT_PATH = Path("stats_ci.json")


def main() -> None:
    print(f"Computing stats for {GFF_SOURCE}")
    stats = compute_gff_stats(GFF_SOURCE)
    if not stats:
        raise SystemExit("compute_gff_stats returned an empty result")

    OUTPUT_PATH.write_text(json.dumps(stats, indent=2))
    print(f"Wrote statistics to {OUTPUT_PATH}")

    schema = json.loads(SCHEMA_PATH.read_text())
    Draft7Validator(schema).validate(stats)
    print(f"Validation succeeded using schema at {SCHEMA_PATH}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise

