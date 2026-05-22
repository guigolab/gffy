"""
Integration test: gffy stats vs Annotrieve production features_statistics.

Fetches the first annotation from the public API, runs ``python -m gffy.cli`` on
``source_file_info.url_path``, and compares the result to ``features_statistics``
stored for that annotation (GFFStats schema in annotrieve embedded documents).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from typing import Any

import pytest

ANNOTRIEVE_API_BASE = os.environ.get(
    "ANNOTRIEVE_API_BASE", "https://genome.crg.es/annotrieve/api/v0"
).rstrip("/")
REQUEST_TIMEOUT = int(os.environ.get("ANNOTRIEVE_API_TIMEOUT", "120"))


def _fetch_first_annotation() -> dict[str, Any]:
    url = f"{ANNOTRIEVE_API_BASE}/annotations?limit=1&offset=5"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "gffy-tests/production-parity"},
    )
    with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
        payload = json.load(resp)
    results = payload.get("results") or []
    if not results:
        pytest.skip("Annotrieve API returned no annotations")
    return results[0]


def _source_url(annotation: dict[str, Any]) -> str | None:
    sfi = annotation.get("source_file_info") or {}
    if isinstance(sfi, dict):
        url_path = sfi.get("url_path")
        if url_path:
            return str(url_path).strip()
    return None


def _run_gffy_cli_stats(source: str) -> dict[str, Any]:
    result = subprocess.run(
        [sys.executable, "-m", "gffy.cli", source],
        capture_output=True,
        text=True,
        timeout=REQUEST_TIMEOUT * 4,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def _assert_stats_equal(production: Any, computed: Any, path: str = "") -> None:
    """Recursively compare GFFStats dicts (production API vs gffy output)."""
    if isinstance(production, dict) and isinstance(computed, dict):
        assert set(production.keys()) == set(computed.keys()), (
            f"key mismatch at {path or '<root>'}: "
            f"prod-only={set(production) - set(computed)}, "
            f"computed-only={set(computed) - set(production)}"
        )
        for key in production:
            child = f"{path}.{key}" if path else key
            _assert_stats_equal(production[key], computed[key], child)
        return

    if isinstance(production, list) and isinstance(computed, list):
        assert len(production) == len(computed), f"list length at {path}"
        for i, (a, b) in enumerate(zip(production, computed)):
            _assert_stats_equal(a, b, f"{path}[{i}]")
        return

    if isinstance(production, float) or isinstance(computed, float):
        assert abs(float(production) - float(computed)) < 1e-9, (
            f"float mismatch at {path}: prod={production!r} computed={computed!r}"
        )
        return

    assert production == computed, (
        f"value mismatch at {path}: prod={production!r} computed={computed!r}"
    )


@pytest.mark.integration
def test_gffy_stats_match_production_annotation():
    """First production annotation: gffy CLI stats == API features_statistics."""
    try:
        annotation = _fetch_first_annotation()
    except urllib.error.URLError as exc:
        pytest.skip(f"Annotrieve API unreachable: {exc}")

    source = _source_url(annotation)
    if not source:
        pytest.skip("First annotation has no source_file_info.url_path")

    production_stats = annotation.get("features_statistics")
    if not production_stats:
        pytest.skip("First annotation has no features_statistics")

    assert "gene_category_stats" in production_stats
    assert "transcript_type_stats" in production_stats

    computed_stats = _run_gffy_cli_stats(source)
    _assert_stats_equal(production_stats, computed_stats)
