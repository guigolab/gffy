"""Shared CLI helpers for gffy and gffy-convert."""

from __future__ import annotations

import argparse
import resource
import sys


def add_usage_argument(parser: argparse.ArgumentParser) -> None:
    """Register ``--usage`` on *parser*."""
    parser.add_argument(
        "--usage",
        action="store_true",
        help="Print user/system time and peak RAM to stderr when the run finishes",
    )


def print_resource_usage(prefix: str, *, stream: object = None) -> None:
    """Print process CPU time and peak RSS for *prefix*."""
    if stream is None:
        stream = sys.stderr
    usage = resource.getrusage(resource.RUSAGE_SELF)
    print(
        f"{prefix} User time: {usage.ru_utime:.2f}s, "
        f"System time: {usage.ru_stime:.2f}s, "
        f"Max RAM: {usage.ru_maxrss / 1024:.1f} MB",
        file=stream,  # type: ignore[arg-type]
    )
