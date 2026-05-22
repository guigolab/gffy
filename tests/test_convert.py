"""Tests for GTF to GFF3 conversion."""

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from gffy import compute_gff_stats, convert_gtf_to_gff3
from gffy.convert import parse_gtf_attributes

_HAS_SORT = shutil.which("sort") is not None


def test_parse_gtf_attributes_quoted_and_unquoted():
    pairs = parse_gtf_attributes('gene_id "G1"; transcript_id T2; tag "a;b"')
    assert pairs == [("gene_id", "G1"), ("transcript_id", "T2"), ("tag", "a;b")]


def test_convert_basic(sample_gtf_path: Path, tmp_path: Path):
    out = tmp_path / "out.gff3"
    summary = convert_gtf_to_gff3(sample_gtf_path, out)
    assert summary["feature_count"] == 12
    assert summary["transcripts_with_cds"] == 1
    assert summary["genes_with_cds"] == 1

    text = out.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert text.startswith("##gff-version 3\n")

    gene1 = next(ln for ln in lines if "\tgene\t" in ln and "ID=G1" in ln)
    assert "biotype=protein_coding" in gene1
    assert "gene_name=GENE1" in gene1
    assert "gene_id" not in gene1

    mrna = next(ln for ln in lines if "\tmRNA\t" in ln)
    assert "ID=T1" in mrna and "Parent=G1" in mrna

    lnc = next(ln for ln in lines if "\tlnc_RNA\t" in ln)
    assert "ID=T2" in lnc and "Parent=G2" in lnc

    pseudo_tx = next(ln for ln in lines if "\ttranscript\t" in ln and "ID=T3" in ln)
    assert "mRNA" not in pseudo_tx.split("\t")[2]

    cds = next(ln for ln in lines if "\tCDS\t" in ln)
    assert "Parent=T1" in cds and "gene_id" not in cds

    exon = next(ln for ln in lines if "\texon\t" in ln and "ID=E1" in ln)
    assert "Parent=T1" in exon

    gene4 = next(ln for ln in lines if "ID=G4" in ln)
    assert "gene_name=weird%3Bname" in gene4


def test_convert_gzip_roundtrip(sample_gtf_path: Path, tmp_path: Path):
    gz_in = tmp_path / "in.gtf.gz"
    with gzip.open(gz_in, "wt", encoding="utf-8") as f:
        f.write(sample_gtf_path.read_text(encoding="utf-8"))

    out = tmp_path / "out.gff3.gz"
    convert_gtf_to_gff3(gz_in, out)
    with gzip.open(out, "rt", encoding="utf-8") as f:
        assert f.readline().strip() == "##gff-version 3"


def test_convert_stats_consumable(sample_gtf_path: Path, tmp_path: Path):
    out = tmp_path / "out.gff3"
    convert_gtf_to_gff3(sample_gtf_path, out)
    stats = compute_gff_stats(str(out))
    coding = stats["gene_category_stats"].get("coding", {})
    assert coding.get("total_count", 0) >= 1


@pytest.mark.skipif(not _HAS_SORT, reason="GNU sort required for low-memory mode")
def test_convert_low_memory_equivalent(sample_gtf_path: Path, tmp_path: Path):
    out_default = tmp_path / "default.gff3"
    out_low = tmp_path / "low.gff3"
    convert_gtf_to_gff3(sample_gtf_path, out_default)
    convert_gtf_to_gff3(
        sample_gtf_path, out_low, force_low_memory=True
    )

    def data_lines(path: Path) -> set[str]:
        return {
            ln
            for ln in path.read_text(encoding="utf-8").splitlines()
            if ln and not ln.startswith("#")
        }

    assert data_lines(out_default) == data_lines(out_low)


def test_convert_cli(sample_gtf_path: Path, tmp_path: Path):
    out = tmp_path / "cli.gff3"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "gffy.convert_cli",
            str(sample_gtf_path),
            "-o",
            str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert out.is_file()
    assert "Mode:" in proc.stderr
    assert proc.stdout == ""


@pytest.mark.skipif(not _HAS_SORT, reason="GNU sort required for low-memory mode")
def test_convert_force_low_memory_mode(sample_gtf_path: Path, tmp_path: Path):
    out = tmp_path / "out.gff3"
    summary = convert_gtf_to_gff3(
        sample_gtf_path, out, force_low_memory=True
    )
    assert "low-memory" in summary["mode"]
