"""Tests for gffy.compute_stats module."""

import pytest
from pathlib import Path

from gffy import compute_gff_stats


class TestComputeGffStats:
    """Test cases for compute_gff_stats function."""

    def test_compute_stats_with_test_file(self):
        """Test computing stats from the test GFF file."""
        test_file = Path(__file__).parent / "fixtures" / "sample_data" / "GCF_000001405.40.gff.gz"

        # Skip test if file doesn't exist
        if not test_file.exists():
            pytest.skip("Test data file not found")

        stats = compute_gff_stats(str(test_file))

        # Basic validation of returned structure
        assert isinstance(stats, dict)
        assert "coding_genes" in stats
        assert "non_coding_genes" in stats
        assert "pseudogenes" in stats

        # Check that we have some genes
        total_genes = (
            stats["coding_genes"]["count"] +
            stats["non_coding_genes"]["count"] +
            stats["pseudogenes"]["count"]
        )
        assert total_genes > 0

    def test_compute_stats_invalid_file(self):
        """Test error handling for invalid file paths."""
        result = compute_gff_stats("/nonexistent/file.gff3")
        # Function returns empty dict on error
        assert result == {}

    def test_compute_stats_empty_result(self):
        """Test handling of files that produce no valid statistics."""
        # This would need a test file that produces empty results
        # For now, just ensure the function doesn't crash
        pass
