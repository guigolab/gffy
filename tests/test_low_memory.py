"""Tests for low-memory (per-seqid + temporary sort) mode."""

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from gffy import compute_gff_stats
from gffy.cache import build_sorted_gff, inspect_source

sort_available = shutil.which("sort") is not None


def _json_keyed(obj: dict) -> str:
    return json.dumps(obj, sort_keys=True)


@pytest.mark.skipif(not sort_available, reason="system sort not available")
class TestLowMemoryMode:
    def test_unsorted_equivalence(self, big_unsorted_gff_path: Path):
        """Per-seqid path must match single-pass on interleaved multi-seqid input."""
        single = compute_gff_stats(str(big_unsorted_gff_path))
        low = compute_gff_stats(
            str(big_unsorted_gff_path),
            force_low_memory=True,
        )
        assert _json_keyed(single) == _json_keyed(low)

    def test_sample_fixture_equivalence(self, sample_gff_path: Path):
        single = compute_gff_stats(str(sample_gff_path))
        low = compute_gff_stats(str(sample_gff_path), force_low_memory=True)
        assert _json_keyed(single) == _json_keyed(low)

    def test_temp_sorted_removed_after_compute(self, sample_gff_path: Path):
        """Sorted file lives only inside a temp dir removed after stats."""
        before = {p.name for p in Path(tempfile.gettempdir()).glob("gffy-*")}
        compute_gff_stats(str(sample_gff_path), force_low_memory=True)
        after = {p.name for p in Path(tempfile.gettempdir()).glob("gffy-*")}
        assert before == after

    def test_build_sorted_direct(self, sample_gff_path: Path, tmp_path: Path):
        info = inspect_source(str(sample_gff_path))
        out = tmp_path / "sorted.gff.gz"
        build_sorted_gff(str(sample_gff_path), out, info)
        assert out.is_file()
        assert out.stat().st_size > 0
