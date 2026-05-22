"""Tests for bulk stats (--from-file)."""

import json
import subprocess
import sys
from io import StringIO
from pathlib import Path

import pytest

from gffy.bulk import compute_bulk_stats, read_source_list


def test_read_source_list_dedup_and_comments(tmp_path: Path):
    path = tmp_path / "sources.txt"
    path.write_text(
        "# header comment\n"
        "/path/a.gff3\n"
        "\n"
        "  /path/b.gff3  \n"
        "/path/a.gff3\n"
        "# another\n"
        "/path/c.gff3\n",
        encoding="utf-8",
    )
    result = read_source_list(str(path))
    assert result.sources == [
        "/path/a.gff3",
        "/path/b.gff3",
        "/path/c.gff3",
    ]
    assert result.duplicates == ["/path/a.gff3"]


def test_read_source_list_stdin(monkeypatch):
    monkeypatch.setattr("sys.stdin", StringIO("/x.gff3\n\n/y.gff3\n/x.gff3\n"))
    result = read_source_list("-")
    assert result.sources == ["/x.gff3", "/y.gff3"]
    assert result.duplicates == ["/x.gff3"]


def test_compute_bulk_stats_success_and_error(
    sample_gff_path: Path,
    big_unsorted_gff_path: Path,
    tmp_path: Path,
):
    missing = tmp_path / "does_not_exist.gff3"
    sources = [str(sample_gff_path), str(big_unsorted_gff_path), str(missing)]
    out = tmp_path / "out.jsonl"

    summary = compute_bulk_stats(sources, out, workers=1)
    assert summary["total"] == 3
    assert summary["succeeded"] == 2
    assert summary["failed"] == 1

    records = [json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]
    assert len(records) == 3
    by_source = {r["source"]: r for r in records}
    assert "error" in by_source[str(missing)]
    assert by_source[str(missing)]["error_category"] == "not_found"
    assert "stats" in by_source[str(sample_gff_path)]
    assert "gene_category_stats" in by_source[str(sample_gff_path)]["stats"]


def test_compute_bulk_stats_workers_parallel(
    sample_gff_path: Path,
    big_unsorted_gff_path: Path,
    tmp_path: Path,
):
    missing = tmp_path / "missing.gff3"
    sources = [str(sample_gff_path), str(big_unsorted_gff_path), str(missing)]
    out = tmp_path / "parallel.jsonl"

    summary = compute_bulk_stats(sources, out, workers=2)
    assert summary["total"] == 3
    records = [json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]
    assert len(records) == 3
    assert len({r["source"] for r in records}) == 3


def test_compute_bulk_stats_fail_fast(tmp_path: Path):
    missing = tmp_path / "missing.gff3"
    out = tmp_path / "failfast.jsonl"

    with pytest.raises(RuntimeError, match="FileNotFoundError"):
        compute_bulk_stats(
            [str(missing)],
            out,
            workers=1,
            continue_on_error=False,
        )

    lines = [ln for ln in out.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert "error" in record


def test_compute_bulk_stats_empty_raises():
    with pytest.raises(ValueError, match="No sources"):
        compute_bulk_stats([], Path("/tmp/x.jsonl"))


def test_cli_bulk_smoke(sample_gff_path: Path, tmp_path: Path):
    list_path = tmp_path / "urls.txt"
    list_path.write_text(f"{sample_gff_path}\n", encoding="utf-8")
    out = tmp_path / "cli.jsonl"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "gffy",
            "--from-file",
            str(list_path),
            "-o",
            str(out),
            "--workers",
            "1",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    records = [json.loads(ln) for ln in out.read_text(encoding="utf-8").splitlines()]
    assert len(records) == 1
    assert "stats" in records[0]


def test_cli_rejects_from_file_without_output(sample_gff_path: Path, tmp_path: Path):
    list_path = tmp_path / "urls.txt"
    list_path.write_text(f"{sample_gff_path}\n", encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, "-m", "gffy", "--from-file", str(list_path)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 2


def test_cli_rejects_both_source_and_from_file(sample_gff_path: Path, tmp_path: Path):
    list_path = tmp_path / "urls.txt"
    list_path.write_text(f"{sample_gff_path}\n", encoding="utf-8")
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "gffy",
            str(sample_gff_path),
            "--from-file",
            str(list_path),
            "-o",
            str(tmp_path / "x.jsonl"),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
