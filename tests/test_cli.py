"""Tests for gffy CLI."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

jsonschema = pytest.importorskip("jsonschema")


class TestCLI:
    def test_cli_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "GFF3 feature statistics" in result.stdout

    def test_python_m_gffy_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "gffy", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "GFF3 feature statistics" in result.stdout

    def test_cli_with_fixture(self, sample_gff_path, schema):
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", str(sample_gff_path)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        stats = json.loads(result.stdout)
        jsonschema.validate(stats, schema)
        assert "gene_category_stats" in stats

    def test_cli_invalid_file(self):
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", "/nonexistent/file.gff3"],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "ERROR (not_found)" in result.stderr

    def test_cli_pretty_output(self, sample_gff_path):
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", str(sample_gff_path), "--pretty"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "\n" in result.stdout
        json.loads(result.stdout)

    def test_cli_output_file(self, sample_gff_path, tmp_path):
        out = tmp_path / "out.json"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                str(sample_gff_path),
                "--output",
                str(out),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert out.exists()
        with open(out, encoding="utf-8") as f:
            stats = json.load(f)
        assert "transcript_type_stats" in stats

    @pytest.mark.skipif(
        __import__("shutil").which("sort") is None,
        reason="system sort not available",
    )
    def test_cli_low_memory(self, sample_gff_path):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                str(sample_gff_path),
                "--low-memory",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "low-memory" in result.stderr
        stats = json.loads(result.stdout)
        assert "gene_category_stats" in stats
