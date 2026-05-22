"""Tests for Annotrieve-compatible CustomAnnotation export."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from gffy.annotrieve_export import (
    build_custom_annotation,
    compute_features_summary_from_lines,
    derive_custom_name,
    format_annotrieve_jsonl_line,
)
from gffy.bulk import compute_bulk_annotrieve_export
from tests.conftest import SAMPLE_GFF3

jsonschema = pytest.importorskip("jsonschema")


@pytest.fixture
def annotrieve_schema() -> dict:
    path = Path(__file__).resolve().parent.parent / "annotrieve-schema.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class TestDeriveCustomName:
    def test_url_uses_path_only(self):
        name = derive_custom_name("https://example.com/path/ann gff3")
        assert name == "path_ann_gff3"
        assert "example.com" not in name

    def test_ftp_url_uses_path_only(self):
        name = derive_custom_name(
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/foo_genomic.gff.gz"
        )
        assert "ncbi.nlm.nih.gov" not in name
        assert name.startswith("genomes_all_GCA")

    def test_override_wins(self):
        assert derive_custom_name("https://x/y", "My experiment") == "My experiment"


class TestFeaturesSummary:
    def test_sample_has_types_and_sources(self):
        lines = [ln for ln in SAMPLE_GFF3.splitlines() if ln and not ln.startswith("#")]
        summary = compute_features_summary_from_lines(lines)
        assert summary["types"]
        assert summary["sources"]
        assert "gene" in summary["types"]
        assert "source" in summary["sources"]
        assert summary["has_cds"] is True
        assert summary["has_exon"] is True


class TestBuildCustomAnnotation:
    def test_small_file_uses_single_pass_stats(self, sample_gff_path, monkeypatch):
        per_seqid_called: list[bool] = []

        def track_per_seqid(path):
            per_seqid_called.append(True)
            from gffy.stats import _compute_per_seqid as real

            return real(path)

        monkeypatch.setattr(
            "gffy.annotrieve_export._compute_per_seqid",
            track_per_seqid,
        )
        build_custom_annotation(str(sample_gff_path))
        assert not per_seqid_called

    def test_force_low_memory_uses_per_seqid(self, sample_gff_path, monkeypatch):
        per_seqid_called: list[bool] = []

        def track_per_seqid(path):
            per_seqid_called.append(True)
            from gffy.stats import _compute_per_seqid as real

            return real(path)

        monkeypatch.setattr(
            "gffy.annotrieve_export._compute_per_seqid",
            track_per_seqid,
        )
        build_custom_annotation(str(sample_gff_path), force_low_memory=True)
        assert per_seqid_called

    def test_e2e_validates_schema(self, sample_gff_path, annotrieve_schema):
        record = build_custom_annotation(str(sample_gff_path))
        jsonschema.validate(record, annotrieve_schema)
        assert record["kind"] == "custom"
        assert len(record["annotation_id"]) == 32
        assert record["annotation_id"] == record["uploaded_md5"]
        assert record["features_statistics"]["gene_category_stats"]

    def test_custom_name_override(self, sample_gff_path):
        record = build_custom_annotation(
            str(sample_gff_path),
            custom_name_override="Lab run 1",
        )
        assert record["custom_name"] == "Lab run 1"

    def test_missing_types_raises(self, tmp_path):
        empty = tmp_path / "empty.gff3"
        empty.write_text("##gff-version 3\n", encoding="utf-8")
        with pytest.raises(ValueError, match="no types or sources"):
            build_custom_annotation(str(empty))


class TestBulkAnnotrieveExport:
    def test_two_ok_one_missing(
        self,
        sample_gff_path,
        big_unsorted_gff_path,
        tmp_path,
    ):
        missing = tmp_path / "missing.gff3"
        out = tmp_path / "import.jsonl"
        summary = compute_bulk_annotrieve_export(
            [str(sample_gff_path), str(big_unsorted_gff_path), str(missing)],
            out,
            workers=1,
        )
        assert summary["succeeded"] == 2
        assert summary["failed"] == 1
        lines = [ln for ln in out.read_text(encoding="utf-8").splitlines() if ln.strip()]
        assert len(lines) == 2
        for line in lines:
            rec = json.loads(line)
            assert rec["kind"] == "custom"

    def test_jsonl_line_format(self, sample_gff_path):
        record = build_custom_annotation(str(sample_gff_path))
        line = format_annotrieve_jsonl_line(record)
        assert "\n" not in line
        parsed = json.loads(line)
        assert parsed["annotation_id"] == record["annotation_id"]


class TestAnnotrieveCLI:
    def test_single_annotrieve_json(self, sample_gff_path, tmp_path, annotrieve_schema):
        out = tmp_path / "fav.json"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                str(sample_gff_path),
                "--annotrieve-json",
                str(out),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert not result.stdout.strip()
        record = json.loads(out.read_text(encoding="utf-8"))
        assert record["kind"] == "custom"
        jsonschema.validate(record, annotrieve_schema)

    def test_single_custom_name(self, sample_gff_path, tmp_path):
        out = tmp_path / "fav.json"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                str(sample_gff_path),
                "--annotrieve-json",
                str(out),
                "--custom-name",
                "My GFF",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        record = json.loads(out.read_text(encoding="utf-8"))
        assert record["custom_name"] == "My GFF"

    def test_bulk_annotrieve_jsonl(self, sample_gff_path, tmp_path):
        list_path = tmp_path / "sources.txt"
        list_path.write_text(f"{sample_gff_path}\n", encoding="utf-8")
        ann_out = tmp_path / "import.jsonl"
        diag_out = tmp_path / "diag.jsonl"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                "--from-file",
                str(list_path),
                "-o",
                str(diag_out),
                "--annotrieve-jsonl",
                str(ann_out),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        assert len([ln for ln in ann_out.read_text().splitlines() if ln.strip()]) == 1
        assert diag_out.exists()

    def test_bulk_bad_source_omitted_from_annotrieve(
        self,
        sample_gff_path,
        tmp_path,
    ):
        missing = tmp_path / "missing.gff3"
        list_path = tmp_path / "sources.txt"
        list_path.write_text(
            f"{sample_gff_path}\n{missing}\n",
            encoding="utf-8",
        )
        ann_out = tmp_path / "import.jsonl"
        diag_out = tmp_path / "diag.jsonl"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gffy.cli",
                "--from-file",
                str(list_path),
                "-o",
                str(diag_out),
                "--annotrieve-jsonl",
                str(ann_out),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert len([ln for ln in ann_out.read_text().splitlines() if ln.strip()]) == 1
        diag_lines = [
            json.loads(ln)
            for ln in diag_out.read_text(encoding="utf-8").splitlines()
            if ln.strip()
        ]
        assert len(diag_lines) == 2
        assert any("error" in r for r in diag_lines)
