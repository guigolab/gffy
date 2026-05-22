"""Tests for error classification and bulk error sidecars."""

import json
from io import StringIO
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
from urllib.error import HTTPError, URLError

import pytest

from gffy.bulk import (
    compute_bulk_stats,
    read_source_list,
    warn_source_list_issues,
)
from gffy.errors import (
    ErrorInfo,
    classify_exception,
    error_record,
    is_fetch_failure,
)


def test_classify_file_not_found():
    info = classify_exception(FileNotFoundError("GFF file not found: x.gff3"))
    assert info.category == "not_found"


def test_classify_http_error():
    exc = HTTPError(
        url="https://example.com/a.gff3",
        code=404,
        msg="Not Found",
        hdrs=None,
        fp=None,
    )
    info = classify_exception(exc)
    assert info.category == "http"
    assert info.http_status == 404


def test_classify_url_error():
    info = classify_exception(URLError("connection refused"))
    assert info.category == "network"
    assert is_fetch_failure(info)


def test_classify_sort_runtime_error():
    info = classify_exception(RuntimeError("sort pipeline failed"))
    assert info.category == "sort"


def test_error_record_includes_category():
    info = ErrorInfo("http", "HTTP 503", http_status=503)
    rec = error_record("https://x/a.gff3", info, mode="single-pass", elapsed_ms=10)
    assert rec["error_category"] == "http"
    assert rec["http_status"] == 503


def test_read_source_list_duplicates(tmp_path: Path, capsys):
    path = tmp_path / "list.txt"
    path.write_text("https://a.gff\nhttps://b.gff\nhttps://a.gff\n", encoding="utf-8")
    result = read_source_list(str(path))
    assert result.sources == ["https://a.gff", "https://b.gff"]
    assert result.duplicates == ["https://a.gff"]
    warn_source_list_issues(result)
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "duplicate" in err


def test_read_source_list_stdin(monkeypatch):
    monkeypatch.setattr("sys.stdin", StringIO("/x.gff3\n/y.gff3\n"))
    result = read_source_list("-")
    assert result.sources == ["/x.gff3", "/y.gff3"]
    assert result.duplicates == []


def test_bulk_not_found_not_in_failed_urls(
    sample_gff_path: Path,
    tmp_path: Path,
):
    missing = tmp_path / "missing.gff3"
    out = tmp_path / "out.jsonl"
    summary = compute_bulk_stats(
        [str(sample_gff_path), str(missing)],
        out,
        workers=1,
    )
    assert summary["failed"] == 1
    failed_urls = Path(str(out) + ".failed_urls.txt")
    assert not failed_urls.exists()

    records = [json.loads(ln) for ln in out.read_text(encoding="utf-8").splitlines()]
    err_rec = next(r for r in records if "error" in r)
    assert err_rec["error_category"] == "not_found"


def test_bulk_failed_urls_sidecar(sample_gff_path: Path, tmp_path: Path):
    bad_url = "https://invalid.example.invalid/annotation.gff3"
    out = tmp_path / "stats.jsonl"

    def fake_run_one(source: str, force_low_memory: bool) -> dict:
        if source == bad_url:
            from gffy.errors import classify_exception, error_record

            exc = URLError("Name or service not known")
            info = classify_exception(exc)
            return error_record(
                source, info, mode="single-pass", elapsed_ms=1, exc=exc
            )
        return {
            "source": source,
            "mode": "single-pass",
            "elapsed_ms": 1,
            "stats": {"gene_category_stats": {}, "transcript_type_stats": {}},
        }

    with patch("gffy.bulk.ProcessPoolExecutor", ThreadPoolExecutor):
        with patch("gffy.bulk._run_one", side_effect=fake_run_one):
            summary = compute_bulk_stats(
                [str(sample_gff_path), bad_url],
                out,
                workers=1,
            )

    sidecar = Path(str(out) + ".failed_urls.txt")
    assert sidecar.is_file()
    assert summary["failed_urls_written"] == str(sidecar)
    urls = sidecar.read_text(encoding="utf-8").strip().splitlines()
    assert bad_url in urls


def test_bulk_custom_failed_urls_out(sample_gff_path: Path, tmp_path: Path):
    bad_url = "https://example.com/missing.gff3"
    out = tmp_path / "out.jsonl"
    custom = tmp_path / "retry.txt"

    def fake_run_one(source: str, force_low_memory: bool) -> dict:
        if source == bad_url:
            from gffy.errors import classify_exception, error_record

            exc = HTTPError(bad_url, 503, "Unavailable", None, None)
            info = classify_exception(exc)
            return error_record(
                source, info, mode="single-pass", elapsed_ms=1, exc=exc
            )
        return {
            "source": source,
            "mode": "single-pass",
            "elapsed_ms": 1,
            "stats": {"gene_category_stats": {}, "transcript_type_stats": {}},
        }

    with patch("gffy.bulk.ProcessPoolExecutor", ThreadPoolExecutor):
        with patch("gffy.bulk._run_one", side_effect=fake_run_one):
            summary = compute_bulk_stats(
                [bad_url],
                out,
                workers=1,
                failed_urls_output=custom,
            )

    assert custom.read_text(encoding="utf-8").strip() == bad_url
    assert summary["failed_urls_written"] == str(custom)
    assert not Path(str(out) + ".failed_urls.txt").exists()


def test_bulk_success_no_failed_urls_file(sample_gff_path: Path, tmp_path: Path):
    out = tmp_path / "ok.jsonl"
    compute_bulk_stats([str(sample_gff_path)], out, workers=1)
    assert not Path(str(out) + ".failed_urls.txt").exists()
