"""
Bulk GFF3 statistics from a list of sources.

Uses a process pool so each file runs in an isolated subprocess with the same
single-pass / low-memory rules as ``compute_gff_stats``.
"""

from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from .annotrieve_export import build_custom_annotation, format_annotrieve_jsonl_line
from .errors import classify_exception, error_record, is_url_source
from .stats import compute_gff_stats, describe_compute_mode

_MAX_WORKERS = 8


@dataclass
class SourceListResult:
    """Parsed bulk input list."""

    sources: list[str] = field(default_factory=list)
    duplicates: list[str] = field(default_factory=list)


def read_source_list(path_or_dash: str) -> SourceListResult:
    """
    Read newline-separated sources from a file or ``-`` (stdin).

    Strips whitespace, skips blank and ``#``-prefixed lines, dedupes while
    preserving first-seen order. Duplicate lines are recorded in ``duplicates``.
    """
    if path_or_dash == "-":
        lines = sys.stdin.readlines()
    else:
        with open(path_or_dash, encoding="utf-8") as fh:
            lines = fh.readlines()

    seen: set[str] = set()
    result = SourceListResult()
    for line in lines:
        entry = line.strip()
        if not entry or entry.startswith("#"):
            continue
        if entry in seen:
            result.duplicates.append(entry)
            continue
        seen.add(entry)
        result.sources.append(entry)
    return result


def warn_source_list_issues(
    result: SourceListResult,
    *,
    prefix: str = "[gffy]",
) -> None:
    """Print duplicate-line warnings to stderr."""
    for dup in result.duplicates:
        print(
            f"{prefix} WARNING: duplicate source (skipped): {dup}",
            file=sys.stderr,
        )
    if result.duplicates:
        print(
            f"{prefix} WARNING: {len(result.duplicates)} duplicate line(s) in input list",
            file=sys.stderr,
        )


def _default_failed_urls_path(output: Path) -> Path:
    return Path(str(output) + ".failed_urls.txt")


def _clamp_workers(workers: int) -> int:
    if workers < 1:
        return 1
    cap = min(os.cpu_count() or 1, _MAX_WORKERS)
    return min(workers, cap)


def _run_one(source: str, force_low_memory: bool) -> dict:
    """Worker entry point (module-level for pickling under spawn)."""
    start = time.monotonic_ns()
    mode_label = "unknown"
    try:
        mode_label, _ = describe_compute_mode(
            source, force_low_memory=force_low_memory
        )
    except Exception:
        pass

    try:
        stats = compute_gff_stats(source, force_low_memory=force_low_memory)
        elapsed_ms = (time.monotonic_ns() - start) // 1_000_000
        return {
            "source": source,
            "mode": mode_label,
            "elapsed_ms": elapsed_ms,
            "stats": stats,
        }
    except Exception as exc:
        elapsed_ms = (time.monotonic_ns() - start) // 1_000_000
        info = classify_exception(exc)
        return error_record(
            source,
            info,
            mode=mode_label,
            elapsed_ms=elapsed_ms,
            exc=exc,
        )


def _run_one_annotrieve(source: str, force_low_memory: bool) -> dict:
    """Worker entry point for Annotrieve CustomAnnotation export."""
    start = time.monotonic_ns()
    try:
        annotation = build_custom_annotation(
            source, force_low_memory=force_low_memory
        )
        elapsed_ms = (time.monotonic_ns() - start) // 1_000_000
        return {
            "ok": True,
            "source": source,
            "elapsed_ms": elapsed_ms,
            "annotation": annotation,
        }
    except Exception as exc:
        elapsed_ms = (time.monotonic_ns() - start) // 1_000_000
        info = classify_exception(exc)
        return {
            "ok": False,
            **error_record(
                source,
                info,
                mode="annotrieve",
                elapsed_ms=elapsed_ms,
                exc=exc,
            ),
        }


def compute_bulk_annotrieve_export(
    sources: Iterable[str],
    output: str | Path,
    *,
    workers: int = 1,
    force_low_memory: bool = False,
    continue_on_error: bool = True,
    progress: Callable[[str, dict], None] | None = None,
) -> dict:
    """
    Build CustomAnnotation records for each source and write import-ready JSONL.

    Only successful records are written to *output*.
    """
    source_list = list(sources)
    if not source_list:
        raise ValueError("No sources to process")

    output = Path(output)
    workers = _clamp_workers(workers)
    total = len(source_list)
    succeeded = 0
    failed = 0

    with open(output, "w", encoding="utf-8", buffering=1) as out_fh:
        executor = ProcessPoolExecutor(max_workers=workers)
        try:
            futures = {
                executor.submit(_run_one_annotrieve, src, force_low_memory): src
                for src in source_list
            }
            for future in as_completed(futures):
                record = future.result()
                if record.get("ok"):
                    succeeded += 1
                    out_fh.write(
                        format_annotrieve_jsonl_line(record["annotation"]) + "\n"
                    )
                    out_fh.flush()
                else:
                    failed += 1
                if progress:
                    progress(record.get("source", ""), record)
                if not record.get("ok") and not continue_on_error:
                    executor.shutdown(wait=False, cancel_futures=True)
                    raise RuntimeError(record.get("error", "export failed"))
        except KeyboardInterrupt:
            executor.shutdown(wait=False, cancel_futures=True)
            raise
        finally:
            executor.shutdown(wait=True)

    return {"total": total, "succeeded": succeeded, "failed": failed}


def compute_bulk_stats(
    sources: Iterable[str],
    output: str | Path,
    *,
    workers: int = 1,
    force_low_memory: bool = False,
    continue_on_error: bool = True,
    failed_urls_output: str | Path | None = None,
    progress: Callable[[str, dict], None] | None = None,
    duplicates_skipped: int = 0,
) -> dict:
    """
    Compute stats for each source concurrently and stream JSONL to *output*.

    Returns summary with ``total``, ``succeeded``, ``failed``,
    ``failed_urls_written``, and ``duplicates_skipped``.
    """
    source_list = list(sources)
    if not source_list:
        raise ValueError("No sources to process")

    output = Path(output)
    workers = _clamp_workers(workers)
    total = len(source_list)
    succeeded = 0
    failed = 0
    failed_url_records: list[str] = []

    with open(output, "w", encoding="utf-8", buffering=1) as out_fh:
        executor = ProcessPoolExecutor(max_workers=workers)
        try:
            futures = {
                executor.submit(_run_one, src, force_low_memory): src
                for src in source_list
            }
            for future in as_completed(futures):
                record = future.result()
                out_fh.write(json.dumps(record, separators=(",", ":")) + "\n")
                out_fh.flush()

                if "error" in record:
                    failed += 1
                    category = record.get("error_category", "unknown")
                    if is_url_source(record["source"]) and category in (
                        "network",
                        "http",
                    ):
                        failed_url_records.append(record["source"])
                    if progress:
                        progress(record["source"], record)
                    if not continue_on_error:
                        executor.shutdown(wait=False, cancel_futures=True)
                        raise RuntimeError(record["error"])
                else:
                    succeeded += 1
                    if progress:
                        progress(record["source"], record)
        except KeyboardInterrupt:
            executor.shutdown(wait=False, cancel_futures=True)
            raise
        finally:
            executor.shutdown(wait=True)

    failed_urls_written: str | None = None
    if failed_url_records:
        urls_path = (
            Path(failed_urls_output)
            if failed_urls_output is not None
            else _default_failed_urls_path(output)
        )
        urls_path.write_text(
            "\n".join(dict.fromkeys(failed_url_records)) + "\n",
            encoding="utf-8",
        )
        failed_urls_written = str(urls_path)
        print(
            f"[gffy] Wrote {len(failed_url_records)} failed URL(s) to {urls_path}",
            file=sys.stderr,
        )

    return {
        "total": total,
        "succeeded": succeeded,
        "failed": failed,
        "failed_urls_written": failed_urls_written,
        "duplicates_skipped": duplicates_skipped,
    }
