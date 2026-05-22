"""
Structured error classification for gffy (stdlib only).

Classification runs only on exception paths; the success path is unchanged.
"""

from __future__ import annotations

import errno
import socket
from dataclasses import dataclass
from typing import Any
from urllib.error import HTTPError, URLError

from ._io import is_url as _is_url


@dataclass(frozen=True)
class ErrorInfo:
    """Classified failure for logging and JSONL output."""

    category: str  # not_found | network | http | sort | parse | unknown
    message: str
    http_status: int | None = None


def classify_exception(exc: BaseException) -> ErrorInfo:
    """Map an exception to a stable category and message."""
    if isinstance(exc, FileNotFoundError):
        return ErrorInfo("not_found", str(exc))

    if isinstance(exc, OSError) and getattr(exc, "errno", None) == errno.ENOENT:
        return ErrorInfo("not_found", str(exc))

    if isinstance(exc, HTTPError):
        status = exc.code if hasattr(exc, "code") else None
        return ErrorInfo(
            "http",
            f"HTTP {status}: {exc.reason}" if status else str(exc),
            http_status=status,
        )

    if isinstance(exc, URLError):
        return ErrorInfo("network", str(exc.reason if exc.reason else exc))

    if isinstance(exc, (TimeoutError, socket.timeout)):
        return ErrorInfo("network", str(exc))

    if isinstance(exc, RuntimeError):
        msg = str(exc).lower()
        if "sort" in msg:
            return ErrorInfo("sort", str(exc))

    return ErrorInfo("unknown", f"{type(exc).__name__}: {exc}")


def is_fetch_failure(info: ErrorInfo) -> bool:
    """True when a URL could not be fetched (network or HTTP error)."""
    return info.category in ("network", "http")


def is_url_source(source: str) -> bool:
    """Return True if *source* is an http(s)/ftp URL."""
    return _is_url(source)


def print_classified_error(
    prefix: str,
    exc: BaseException,
    *,
    stream: Any = None,
) -> None:
    """Print a classified error message to *stream* (default stderr)."""
    import sys

    if stream is None:
        stream = sys.stderr
    info = classify_exception(exc)
    print(f"{prefix} ERROR ({info.category}): {info.message}", file=stream)
    if info.http_status is not None:
        print(f"{prefix} HTTP status: {info.http_status}", file=stream)


def error_record(
    source: str,
    info: ErrorInfo,
    *,
    mode: str,
    elapsed_ms: int,
    exc: BaseException | None = None,
) -> dict[str, Any]:
    """Build a JSONL failure record."""
    if exc is not None:
        error_text = f"{type(exc).__name__}: {exc}"
    else:
        error_text = info.message

    record: dict[str, Any] = {
        "source": source,
        "mode": mode,
        "elapsed_ms": elapsed_ms,
        "error": error_text,
        "error_category": info.category,
    }
    if info.http_status is not None:
        record["http_status"] = info.http_status
    return record
