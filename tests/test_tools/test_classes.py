"""Tests for gffy.tools.classes module."""

from array import array
from gffy.tools.classes import Root, Transcript, OrphanFeature


class TestRoot:
    """Test cases for Root class."""

    def test_root_initialization(self):
        """Test Root initialization."""
        root = Root("gene", "protein_coding", 1000)

        assert root.feature_type == "gene"
        assert root.biotype == "protein_coding"
        assert root.length == 1000
        assert root.has_cds is False
        assert root.has_exon is False

    def test_root_with_none_biotype(self):
        """Test Root with None biotype."""
        root = Root("gene", None, 500)

        assert root.biotype is None
        assert root.length == 500


class TestTranscript:
    """Test cases for Transcript class."""

    def test_transcript_initialization(self):
        """Test Transcript initialization."""
        transcript = Transcript("gene123", "mRNA")

        assert transcript.gene_id == "gene123"
        assert transcript.type == "mRNA"
        assert isinstance(transcript.exons_flat, array)
        assert transcript.exon_len_sum == 0
        assert transcript.cds_total_len == 0
        assert transcript.cds_segments == 0
        assert isinstance(transcript.cds_lens, array)

    def test_transcript_with_none_type(self):
        """Test Transcript with None type."""
        transcript = Transcript("gene456", None)

        assert transcript.type is None


class TestOrphanFeature:
    """Test cases for OrphanFeature class."""

    def test_orphan_feature_initialization(self):
        """Test OrphanFeature initialization."""
        orphan = OrphanFeature("exon123", "exon", 100, 200, ["transcript456"], "protein_coding")

        assert orphan.feature_id == "exon123"
        assert orphan.feature_type == "exon"
        assert orphan.start == 100
        assert orphan.end == 200
        assert orphan.parent_ids == ["transcript456"]
        assert orphan.biotype == "protein_coding"

    def test_orphan_feature_with_none_biotype(self):
        """Test OrphanFeature with None biotype."""
        orphan = OrphanFeature("cds789", "CDS", 50, 150, ["transcript123", "gene456"], None)

        assert orphan.biotype is None
        assert len(orphan.parent_ids) == 2
