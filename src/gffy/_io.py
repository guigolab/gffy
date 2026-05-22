"""Private I/O helpers shared across gffy modules."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator, TextIO

_GZ_SUFFIXES = (".gz", ".gzip")
_URL_SCHEMES = ("http://", "https://", "ftp://")


def is_gzipped_path(path: str | Path) -> bool:
    """Return True if *path* has a gzip filename suffix."""
    lower = str(path).lower()
    return lower.endswith(_GZ_SUFFIXES)


def is_url(source: str) -> bool:
    """Return True if *source* is an http(s)/ftp URL."""
    return source.startswith(_URL_SCHEMES)


def open_text_auto(path: str | Path, mode: str = "r") -> TextIO:
    """Open a local plain or gzipped text file in text mode."""
    path = str(path)
    if is_gzipped_path(path):
        gz_mode = mode if "t" in mode else f"{mode}t"
        return gzip.open(path, gz_mode, encoding="utf-8", errors="replace")  # type: ignore[return-value]
    plain_mode = mode.replace("t", "") if "t" in mode else mode
    return open(path, plain_mode, encoding="utf-8", errors="replace")


def iter_text_lines(
    path: str | Path,
    *,
    skip_comments: bool = True,
) -> Iterator[str]:
    """Yield non-empty data lines from a local plain or gzipped text file."""
    with open_text_auto(path) as fh:
        for line in fh:
            if skip_comments and (line.startswith("#") or not line.strip()):
                continue
            yield line.rstrip("\n\r")
