"""Shared pytest fixtures for gffy tests."""

import json
from pathlib import Path

import pytest

SAMPLE_GFF3 = """##gff-version 3
chr1\tsource\tgene\t1000\t5000\t.\t+\t.\tID=gene1;biotype=protein_coding
chr1\tsource\tmRNA\t1100\t4900\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding
chr1\tsource\texon\t1100\t1200\t.\t+\t.\tParent=tx1
chr1\tsource\texon\t2000\t2100\t.\t+\t.\tParent=tx1
chr1\tsource\tCDS\t1100\t1200\t.\t+\t.\tParent=tx1
chr1\tsource\tCDS\t2000\t2100\t.\t+\t.\tParent=tx1
chr1\tsource\tpseudogene\t6000\t7000\t.\t+\t.\tID=gene2;biotype=pseudogene
chr1\tsource\ttranscript\t6100\t6900\t.\t+\t.\tID=tx2;Parent=gene2
chr1\tsource\texon\t6100\t6200\t.\t+\t.\tParent=tx2
chr1\tsource\tgene\t8000\t9000\t.\t+\t.\tID=gene3;biotype=lncRNA
chr1\tsource\tlncRNA\t8100\t8900\t.\t+\t.\tID=tx3;Parent=gene3;biotype=lncRNA
chr1\tsource\texon\t8100\t8200\t.\t+\t.\tParent=tx3
chr1\tsource\tgene\t10000\t11000\t.\t+\t.\tID=gene4
chr1\tsource\tmRNA\t10100\t10900\t.\t+\t.\tID=tx4;Parent=gene4
chr1\tsource\texon\t10100\t10200\t.\t+\t.\tParent=tx4
"""

# Records interleaved across seqids (requires sort for per-seqid mode equivalence)
UNSORTED_GFF3 = """##gff-version 3
chr2\tsource\tgene\t1000\t2000\t.\t+\t.\tID=geneB;biotype=protein_coding
chr1\tsource\tgene\t1000\t5000\t.\t+\t.\tID=geneA;biotype=protein_coding
chr3\tsource\tgene\t1000\t3000\t.\t+\t.\tID=geneC;biotype=lncRNA
chr1\tsource\tmRNA\t1100\t4900\t.\t+\t.\tID=txA;Parent=geneA;biotype=protein_coding
chr2\tsource\tmRNA\t1100\t1900\t.\t+\t.\tID=txB;Parent=geneB;biotype=protein_coding
chr1\tsource\texon\t1100\t1200\t.\t+\t.\tParent=txA
chr3\tsource\tlncRNA\t1100\t2900\t.\t+\t.\tID=txC;Parent=geneC;biotype=lncRNA
chr2\tsource\texon\t1100\t1200\t.\t+\t.\tParent=txB
chr1\tsource\tCDS\t1100\t1200\t.\t+\t.\tParent=txA
chr3\tsource\texon\t1100\t1200\t.\t+\t.\tParent=txC
chr2\tsource\tCDS\t1100\t1200\t.\t+\t.\tParent=txB
"""

SAMPLE_GTF = """##description test
chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_version "1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1\tEnsembl\ttranscript\t1100\t4900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; transcript_name "TX1"; transcript_biotype "protein_coding";
chr1\tEnsembl\texon\t1100\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_id "E1"; exon_number "1";
chr1\tEnsembl\tCDS\t1100\t1200\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "1"; ccds_id "CCDS1";
chr1\tEnsembl\tCDS\t2000\t2100\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "2";
chr1\tEnsembl\tgene\t6000\t7000\t.\t+\t.\tgene_id "G2"; gene_name "GENE2"; gene_biotype "lncRNA";
chr1\tEnsembl\tlnc_RNA\t6100\t6900\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; transcript_biotype "lncRNA";
chr1\tEnsembl\texon\t6100\t6200\t.\t+\t.\tgene_id "G2"; transcript_id "T2"; exon_id "E2";
chr1\tEnsembl\tgene\t8000\t9000\t.\t+\t.\tgene_id "G3"; gene_biotype "pseudogene";
chr1\tEnsembl\ttranscript\t8100\t8900\t.\t+\t.\tgene_id "G3"; transcript_id "T3";
chr1\tEnsembl\texon\t8100\t8200\t.\t+\t.\tgene_id "G3"; transcript_id "T3";
chr1\tEnsembl\tgene\t10000\t11000\t.\t+\t.\tgene_id "G4"; gene_name "weird;name";
"""


@pytest.fixture
def sample_gff_path(tmp_path: Path) -> Path:
    """Minimal GFF3 with coding, pseudogene, non-coding, and missing-biotype genes."""
    path = tmp_path / "sample.gff3"
    path.write_text(SAMPLE_GFF3, encoding="utf-8")
    return path


@pytest.fixture
def sample_gtf_path(tmp_path: Path) -> Path:
    """GTF with coding, lncRNA, and pseudogene (no CDS on T3)."""
    path = tmp_path / "sample.gtf"
    path.write_text(SAMPLE_GTF, encoding="utf-8")
    return path


@pytest.fixture
def big_unsorted_gff_path(tmp_path: Path) -> Path:
    """GFF3 with features interleaved across multiple seqids."""
    path = tmp_path / "unsorted.gff3"
    path.write_text(UNSORTED_GFF3, encoding="utf-8")
    return path


@pytest.fixture
def schema() -> dict:
    """Load schema.json from package root."""
    schema_path = Path(__file__).resolve().parent.parent / "schema.json"
    with open(schema_path, encoding="utf-8") as f:
        return json.load(f)
