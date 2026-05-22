"""Tests for gffy.stats module."""

import json

import pytest

jsonschema = pytest.importorskip("jsonschema")

from gffy import compute_gff_stats, compute_gff_stats_from_lines
from gffy.stats import MISSING_BIOTYPE


class TestComputeGffStats:
    def test_top_level_keys(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        assert "gene_category_stats" in stats
        assert "transcript_type_stats" in stats

    def test_validates_against_schema(self, sample_gff_path, schema):
        stats = compute_gff_stats(str(sample_gff_path))
        jsonschema.validate(stats, schema)

    def test_gene_categories(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        cats = stats["gene_category_stats"]
        assert cats["coding"]["total_count"] == 1
        assert cats["pseudogene"]["total_count"] == 1
        assert cats["non_coding"]["total_count"] == 2

    def test_missing_biotype_fallback(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        nc_biotypes = stats["gene_category_stats"]["non_coding"]["biotype_counts"]
        assert MISSING_BIOTYPE in nc_biotypes

    def test_transcript_type_ordering(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        types = list(stats["transcript_type_stats"].keys())
        counts = [stats["transcript_type_stats"][t]["total_count"] for t in types]
        assert counts == sorted(counts, reverse=True)

    def test_exon_concatenated_length_when_multi_exon(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        mrna = stats["transcript_type_stats"]["mRNA"]
        assert "concatenated_length" in mrna["exon_stats"]
        assert mrna["exon_stats"]["total_count"] > mrna["total_count"]

    def test_single_exon_no_concatenated_length(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        # tx4 has one exon; lncRNA and transcript types may differ
        for ttype, tstats in stats["transcript_type_stats"].items():
            if tstats["exon_stats"]["total_count"] <= tstats["total_count"]:
                assert "concatenated_length" not in tstats["exon_stats"]

    def test_cds_stats_only_when_cds_present(self, sample_gff_path):
        stats = compute_gff_stats(str(sample_gff_path))
        assert "cds_stats" in stats["transcript_type_stats"]["mRNA"]
        assert "cds_stats" not in stats["transcript_type_stats"].get("lncRNA", {})

    def test_from_lines_unordered(self, sample_gff_path):
        lines = sample_gff_path.read_text(encoding="utf-8").splitlines()
        # Reverse line order to simulate unordered input
        data_lines = [ln for ln in lines if not ln.startswith("#")]
        reversed_lines = [ln for ln in lines if ln.startswith("#")] + list(reversed(data_lines))
        stats = compute_gff_stats_from_lines(reversed_lines)
        assert stats["gene_category_stats"]["coding"]["total_count"] == 1

    def test_invalid_file_raises(self):
        with pytest.raises(FileNotFoundError):
            compute_gff_stats("/nonexistent/file.gff3")
