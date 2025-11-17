"""Tests for gffy CLI."""

import json
import subprocess
import sys
from pathlib import Path
import pytest


class TestCLI:
    """Test cases for the gffy command-line interface."""

    def test_cli_help(self):
        """Test that CLI shows help message."""
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Compute comprehensive statistics" in result.stdout

    def test_cli_with_test_file(self):
        """Test CLI with the test GFF file."""
        test_file = Path(__file__).parent / "fixtures" / "sample_data" / "GCF_000001405.40.gff.gz"

        if not test_file.exists():
            pytest.skip("Test data file not found")

        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", str(test_file)],
            capture_output=True,
            text=True
        )

        # Should succeed and produce JSON output
        assert result.returncode == 0
        assert result.stderr.strip()  # Should have progress messages on stderr

        # Parse JSON output from stdout
        try:
            stats = json.loads(result.stdout)
            assert isinstance(stats, dict)
            assert "coding_genes" in stats
        except json.JSONDecodeError:
            pytest.fail("CLI did not produce valid JSON output")

    def test_cli_invalid_file(self):
        """Test CLI error handling for invalid files."""
        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", "/nonexistent/file.gff3"],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0
        assert "Error:" in result.stderr

    def test_cli_pretty_output(self):
        """Test CLI pretty-printed JSON output."""
        test_file = Path(__file__).parent / "fixtures" / "sample_data" / "GCF_000001405.40.gff.gz"

        if not test_file.exists():
            pytest.skip("Test data file not found")

        result = subprocess.run(
            [sys.executable, "-m", "gffy.cli", str(test_file), "--pretty"],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0
        # Pretty output should have indentation
        assert "  " in result.stdout or "\n" in result.stdout
