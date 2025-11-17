"""Tests for gffy.tools.helpers module."""

from gffy.tools.helpers import (
    init_orphan_feature,
    categorize_roots,
    init_category_totals,
    parse_gff_line_fast,
    process_feature,
    resolve_orphans,
)
from gffy.tools.classes import Root, OrphanFeature


class TestInitOrphanFeature:
    """Test cases for init_orphan_feature function."""

    def test_init_orphan_feature_basic(self):
        """Test basic orphan feature initialization."""
        orphan = init_orphan_feature(
            "exon123", "exon", 100, 200, ["transcript456"], "protein_coding"
        )

        assert orphan.feature_id == "exon123"
        assert orphan.feature_type == "exon"
        assert orphan.start == 100
        assert orphan.end == 200
        assert orphan.parent_ids == ["transcript456"]
        assert orphan.biotype == "protein_coding"

    def test_init_orphan_feature_none_biotype(self):
        """Test orphan feature with None biotype."""
        orphan = init_orphan_feature("cds789", "CDS", 50, 150, ["transcript123"], None)

        assert orphan.biotype is None


class TestCategorizeRoots:
    """Test cases for categorize_roots function."""

    def test_categorize_pseudogene(self):
        """Test categorization of pseudogenes."""
        roots = {
            "gene1": Root("pseudogene", None, 1000)
        }
        categories = categorize_roots(roots)

        assert categories["gene1"] == "pseudogene"

    def test_categorize_coding_gene(self):
        """Test categorization of coding genes."""
        roots = {
            "gene1": Root("gene", "protein_coding", 1000)
        }
        roots["gene1"].has_cds = True

        categories = categorize_roots(roots)

        assert categories["gene1"] == "coding"

    def test_categorize_non_coding_gene(self):
        """Test categorization of non-coding genes."""
        roots = {
            "gene1": Root("gene", "lncRNA", 1000)
        }
        roots["gene1"].has_exon = True

        categories = categorize_roots(roots)

        assert categories["gene1"] == "non_coding"

    def test_categorize_uncategorized_gene(self):
        """Test categorization of genes that don't fit other categories."""
        roots = {
            "gene1": Root("gene", None, 1000)
        }
        categories = categorize_roots(roots)

        assert categories["gene1"] is None


class TestInitCategoryTotals:
    """Test cases for init_category_totals function."""

    def test_init_category_totals_structure(self):
        """Test that init_category_totals returns correct structure."""
        totals = init_category_totals()

        expected_categories = ["coding", "non_coding", "pseudogene"]
        expected_keys = [
            "exons", "introns", "cds", "exon_len_sum",
            "intron_len_sum", "cds_len_sum", "cds_transcripts"
        ]

        for category in expected_categories:
            assert category in totals
            for key in expected_keys:
                assert key in totals[category]
                assert totals[category][key] == 0


class TestParseGffLineFast:
    """Test cases for parse_gff_line_fast function."""

    def test_parse_basic_gff_line(self):
        """Test parsing a basic GFF line."""
        line = [
            "chr1", "source", "gene", "100", "200", ".", "+", ".",
            "ID=gene1;biotype=protein_coding"
        ]

        feature_type, start, end, parent_ids, biotype, feature_id = parse_gff_line_fast(line)

        assert feature_type == "gene"
        assert start == 100
        assert end == 200
        assert parent_ids == []
        assert biotype == "protein_coding"
        assert feature_id == "gene1"

    def test_parse_gff_line_with_parent(self):
        """Test parsing GFF line with parent relationships."""
        line = [
            "chr1", "source", "exon", "150", "180", ".", "+", ".",
            "ID=exon1;Parent=transcript1,gene1;biotype=protein_coding"
        ]

        feature_type, start, end, parent_ids, biotype, feature_id = parse_gff_line_fast(line)

        assert feature_type == "exon"
        assert parent_ids == ["transcript1", "gene1"]
        assert feature_id == "exon1"


class TestProcessFeature:
    """Test cases for process_feature function."""

    def test_process_root_feature(self):
        """Test processing a root feature (gene)."""
        roots = {}
        id_to_root = {}
        transcripts = {}

        result = process_feature(
            "gene1", "gene", 100, 200, 101, [], "protein_coding",
            roots, id_to_root, transcripts
        )

        assert result is True
        assert "gene1" in roots
        assert roots["gene1"].feature_type == "gene"
        assert id_to_root["gene1"] == "gene1"

    def test_process_transcript_feature(self):
        """Test processing a transcript feature."""
        roots = {"gene1": Root("gene", "protein_coding", 1000)}
        id_to_root = {"gene1": "gene1"}
        transcripts = {}

        result = process_feature(
            "transcript1", "mRNA", 100, 200, 101, ["gene1"], None,
            roots, id_to_root, transcripts
        )

        assert result is True
        assert "transcript1" in transcripts
        assert transcripts["transcript1"].gene_id == "gene1"
        assert transcripts["transcript1"].type == "mRNA"


class TestResolveOrphans:
    """Test cases for resolve_orphans function."""

    def test_resolve_orphans_basic(self):
        """Test basic orphan resolution."""
        roots = {"gene1": Root("gene", "protein_coding", 1000)}
        id_to_root = {"gene1": "gene1"}
        transcripts = {}

        orphans = [
            OrphanFeature("exon1", "exon", 100, 150, ["transcript1"], None)
        ]

        # First add transcript1 to id_to_root
        id_to_root["transcript1"] = "gene1"

        still_orphaned = resolve_orphans(orphans, roots, id_to_root, transcripts)

        # Should process successfully
        assert len(still_orphaned) == 0
        assert "transcript1" in transcripts
