# gffy

[![CI](https://github.com/guigolab/gffy/actions/workflows/ci.yml/badge.svg)](https://github.com/guigolab/gffy/actions/workflows/ci.yml)
[![Publish to PyPI](https://github.com/guigolab/gffy/actions/workflows/publish.yml/badge.svg)](https://github.com/guigolab/gffy/actions/workflows/publish.yml)
[![PyPI version](https://img.shields.io/pypi/v/gffy.svg)](https://pypi.org/project/gffy/)
[![PyPI downloads](https://img.shields.io/pypi/dm/gffy.svg)](https://pypi.org/project/gffy/)
[![Python versions](https://img.shields.io/pypi/pyversions/gffy.svg)](https://pypi.org/project/gffy/)
[![PyPI format](https://img.shields.io/pypi/format/gffy.svg)](https://pypi.org/project/gffy/)
[![License](https://img.shields.io/pypi/l/gffy.svg)](https://pypi.org/project/gffy/)

Fast, low-RAM GFF3 feature statistics aligned with the **annotrieve** `GFFStats` schema (`gene_category_stats` and `transcript_type_stats`), plus a **GTF → GFF3** converter (`gffy-convert`).

Use **gffy** for statistics JSON, **gffy-convert** for format conversion, or import either workflow in Python pipelines and services.

## Installation

### From PyPI

```bash
pip install gffy
```

### From source (development)

```bash
git clone https://github.com/guigolab/gffy.git
cd gffy
pip install -e ".[dev]"
```

**Runtime:** Python 3.9+, stdlib only. **Low-memory mode** additionally requires GNU `sort` on `PATH`.

## Command-line usage

After installation, the `gffy` console script is on your PATH. You can also run:

```bash
python -m gffy
```

### Options

| Argument / flag | Description |
|-----------------|-------------|
| `gff_source` | Local path or URL to a GFF3 file (plain or `.gz`); mutually exclusive with `--from-file` |
| `--from-file` | Newline-separated list of paths/URLs (`-` = stdin); writes **JSONL** (requires `-o`) |
| `-o`, `--output` | Write JSON (single file) or JSONL (bulk) to this path |
| `--pretty` | Pretty-print JSON (single-source only; not valid with `--from-file`) |
| `--low-memory` | Force sort-to-temp + per-seqid mode (applies to every source in bulk) |
| `--workers N` | Parallel worker processes for `--from-file` (default `1`, max `min(cpu_count, 8)`) |
| `--fail-fast` | Stop bulk run on first failing source (default: record error and continue) |
| `--failed-urls-out` | Override path for unfetchable URL list (default: `<output>.failed_urls.txt`) |
| `--usage` | Print user/system time and peak RAM to stderr at end of run |
| `--annotrieve-json` | Write one **Annotrieve Favorites** import JSON file (single source) |
| `--custom-name` | Override `custom_name` for `--annotrieve-json` (single source only) |
| `--annotrieve-jsonl` | Write import-ready JSONL for `--from-file` bulk (one `CustomAnnotation` per line) |

### Annotrieve Favorites import (offline)

Compute metadata for [Annotrieve Favorites](https://genome.crg.es/annotrieve/annotations/favorites/) **without uploading GFF** to the server (avoids the daily GFF upload limit). Output matches the **Import from JSON** / **Import from JSONL** drawer: `kind`, `annotation_id`, `uploaded_md5` (MD5 of sorted uncompressed GFF bytes), `features_summary`, and `features_statistics`.

```bash
# Single annotation → JSON file
gffy annotation.gff3 --annotrieve-json my-annotation.json

# Optional display name
gffy annotation.gff3 --annotrieve-json out.json --custom-name "My experiment"

# Stats and Annotrieve export together
gffy annotation.gff3 -o stats.json --annotrieve-json favorites.json

# Bulk: diagnostic JSONL (-o) plus import-ready library file
gffy --from-file urls.txt -o run.jsonl --annotrieve-jsonl favorites-import.jsonl
```

In the UI: **Favorites → Add custom → Import from JSON** or **Import from JSONL**. Re-importing the same sorted content updates the entry by MD5. Schema: [`annotrieve-schema.json`](annotrieve-schema.json) (shipped in the PyPI package).

### Examples

```bash
# JSON to stdout (progress and RAM on stderr)
gffy annotation.gff3

# Save formatted JSON
gffy annotation.gff3.gz --output stats.json --pretty

# Remote URL
gffy https://example.com/annotation.gff3.gz -o stats.json

# Large files: sort by seqid into a temp file, stream per chromosome, then delete temp
gffy large.gff.gz --low-memory
```

**Streams:** statistics JSON goes to **stdout** (or `--output`). Mode label and progress go to **stderr**; add `--usage` for timing and peak RAM at the end.

### Bulk stats from a list

Provide a text file with one path or URL per line (plain or `.gz`). Blank lines and lines starting with `#` are ignored. Duplicate lines are skipped (first occurrence kept); duplicates trigger **warnings** on stderr.

```bash
gffy --from-file urls.txt -o stats.jsonl
gffy --from-file urls.txt -o stats.jsonl --workers 4 --low-memory
cat urls.txt | gffy --from-file - -o stats.jsonl
gffy --from-file urls.txt -o stats.jsonl --fail-fast
```

Example `urls.txt`:

```text
# Ensembl release
https://example.com/annotation1.gff3.gz
/path/to/local/annotation2.gff3
```

**JSONL output** — one compact JSON object per line (not a JSON array):

Success:

```json
{"source": "https://example.com/a.gff3.gz", "mode": "single-pass", "elapsed_ms": 1234, "stats": {"gene_category_stats": {}, "transcript_type_stats": {}}}
```

Failure (when not using `--fail-fast`):

```json
{"source": "/missing.gff3", "mode": "single-pass", "elapsed_ms": 12, "error": "FileNotFoundError: GFF file not found: /missing.gff3", "error_category": "not_found"}
```

HTTP / connection failures on URLs also set `error_category` (`http` or `network`) and optional `http_status`.

Each `stats` object matches the annotrieve `GFFStats` schema ([`schema.json`](schema.json)).

### Error handling

Errors are classified without extra I/O on the success path (classification runs only when an exception occurs).

| `error_category` | Typical cause |
|------------------|----------------|
| `not_found` | Local path does not exist |
| `network` | Connection refused, timeout, DNS failure |
| `http` | HTTP 4xx/5xx from server (`http_status` field set) |
| `sort` | GNU `sort` pipeline failed in low-memory mode |
| `unknown` | Other failures |

**Bulk (`--from-file`):** failures are **non-blocking** by default — each source gets one JSONL line (success or error). Duplicate entries in the input list print `[gffy] WARNING: duplicate source (skipped): ...` on stderr.

**Failed URL retry list:** when a URL fails with `network` or `http`, it is appended to `<output>.failed_urls.txt` (e.g. `stats.jsonl` → `stats.jsonl.failed_urls.txt`). Use that file to retry only the unfetchable URLs. Local `not_found` paths are not included. Override with `--failed-urls-out PATH`.

**Single-source `gffy` and `gffy-convert`:** failures print `[gffy] ERROR (<category>): ...` on stderr and exit with code 1 (no partial output).

**Memory:** bulk mode runs each source in a **separate process** (`--workers`). Every worker uses the same single-pass / low-memory rules as a normal `gffy` run. Approximate peak RAM scales with `workers × per-file usage` (~500 MB target per file). Default `--workers 1` keeps memory predictable; increase for throughput when you have spare RAM.

**Python API:**

```python
from gffy import compute_bulk_stats, read_source_list, warn_source_list_issues

parsed = read_source_list("urls.txt")
warn_source_list_issues(parsed)
summary = compute_bulk_stats(
    parsed.sources,
    "stats.jsonl",
    workers=2,
    force_low_memory=False,
    continue_on_error=True,
)
print(summary)
# {"total": N, "succeeded": K, "failed": M, "failed_urls_written": "...", "duplicates_skipped": D}
```

## GTF → GFF3 conversion

Use **`gffy-convert`** to translate GTF annotations into GFF3 with `ID` / `Parent` attributes (separate from the stats tool):

```bash
gffy-convert annotation.gtf -o annotation.gff3
gffy-convert annotation.gtf.gz -o annotation.gff3.gz
gffy-convert https://example.com/annotation.gtf.gz -o local.gff3 --low-memory
gffy-convert annotation.gtf -o out.gff3 --no-sort
```

You can also run:

```bash
python -m gffy.convert_cli annotation.gtf -o annotation.gff3
```

### `gffy-convert` options

| Argument / flag | Description |
|-----------------|-------------|
| `gtf_source` | Local path or URL to a GTF file (plain or `.gz`) |
| `-o`, `--output` | **Required** output GFF3 path (use `.gz` for gzip) |
| `--low-memory` | Force sort-by-seqid into a temp file, then convert |
| `--no-sort` | Keep input row order (skip seqid sort) |
| `--usage` | Print user/system time and peak RAM to stderr at end of run |

**Streams:** GFF3 goes to `--output`. Mode label and feature counts go to **stderr**; add `--usage` for timing and peak RAM at the end.

### Conversion rules

| Input (GTF) | Output (GFF3) |
|---------------|---------------|
| `gene_id "G1"` on a root feature | `ID=G1` |
| `transcript_id "T1"` + `gene_id "G1"` on a transcript row | `ID=T1;Parent=G1` |
| `exon` / `CDS` / UTR rows with `transcript_id` | `Parent=T1` (and `ID=` when `exon_id` is present) |
| Feature type `transcript` **with** CDS children | Type rewritten to **`mRNA`** |
| Feature type `transcript` **without** CDS | Type stays **`transcript`** |
| Root feature whose transcripts have CDS | `biotype=protein_coding` added if no `biotype` already set |

GTF-only keys (`gene_id`, `transcript_id`, `exon_id`) are not copied verbatim; other attributes (e.g. `gene_name`, `transcript_biotype`, `ccds_id`) are percent-encoded in GFF3 form.

Example:

```text
# GTF
chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_biotype "protein_coding";
chr1\tEnsembl\ttranscript\t1100\t4900\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
chr1\tEnsembl\tCDS\t1100\t1200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";

# GFF3 (abbreviated column 9)
chr1\t...\tgene\t...\tID=G1;biotype=protein_coding;gene_biotype=protein_coding
chr1\t...\tmRNA\t...\tID=T1;Parent=G1
chr1\t...\tCDS\t...\tParent=T1
```

### Output row order

`gffy-convert` does **not** reorder features during conversion itself. The GFF3 file is written in the **same line order** as the input stream used for pass 2 (which matches pass 1). Whether that stream is the original GTF or a sorted copy depends on file size and flags.

**When order matches the original GTF**

- **Local file**, below the size thresholds (gz ≤ 100 MB, plain ≤ 1 GB — same as `gffy` stats), and **without** `--low-memory`: the tool reads your file twice in place and emits GFF3 rows in **identical order** to the input (aside from the added `##gff-version 3` header).
- **Any input** with **`--no-sort`**: sorting is disabled even for large files or with `--low-memory`. Row order follows the staged/downloaded file as-is.

**When rows are reordered (seqid sort)**

Sorting runs only when **both** are true:

1. Sorting is enabled (default; turned off with `--no-sort`), and  
2. **`--low-memory`** is set, **or** the source is considered **big** (gz &gt; 100 MB or plain &gt; 1 GB by `Content-Length` / file size).

In that case the GTF is staged and run through GNU `sort` with **`-k1,1`** (first column = seqid only). Effects:

- All rows for a given seqid are grouped together.
- **Across seqids**, order becomes lexicographic by seqid (e.g. `chr1`, `chr10`, `chr2` under `LC_ALL=C`).
- **Within a seqid**, order is preserved relative to the input because sort uses **`-s` (stable)** — gene / transcript / exon / CDS lines on the same chromosome keep their previous relative order.

The converter never sorts by start coordinate, feature type, or hierarchy; only seqid grouping is applied.

**URLs**

Remote sources are always downloaded once into a temporary file (so two passes do not hit the network twice). That temp file is sorted only under the rules above; otherwise output order matches the download order.

**Pass 1 vs pass 2**

Pass 1 scans the file to find which transcripts have CDS (for `mRNA` / `biotype=protein_coding` rules). Pass 2 writes GFF3 lines in the **same order** as the same input stream. No second shuffle happens between passes.

**Quick reference**

| Scenario | Output line order |
|----------|-------------------|
| Small local file, default flags | Same as input GTF |
| Large local file, default flags | Sorted by seqid (stable within seqid) |
| Any size + `--low-memory` | Sorted by seqid (unless `--no-sort`) |
| Any size + `--no-sort` | Same as input (or download order for URLs) |
| Small URL, default flags | Same as download order (no sort) |
| Large URL, default flags | Sorted by seqid |

Use **`--no-sort`** when downstream tools require the exact row order of the source GTF (e.g. diffing against the original). Use the default (or **`--low-memory`** on smaller files) when you want seqids grouped for streaming or parity with `gffy` stats low-memory mode.

### Conversion memory modes

| Mode | When | Behavior |
|------|------|----------|
| **single-pass** | Small files (same size thresholds as stats) | Two passes over the source path; **input row order preserved** |
| **low-memory** | Large files, `--low-memory`, or URL staging + sort | Download/stage once, **sort by seqid** into a **temporary** file, convert, delete temp |

Nothing is cached permanently on disk.

### Python API (conversion)

```python
from gffy import convert_gtf_to_gff3

summary = convert_gtf_to_gff3(
    "/path/to/annotation.gtf.gz",
    "out.gff3.gz",
    force_low_memory=False,
    sort=True,
)
print(summary["feature_count"], summary["genes_with_cds"])
```

| Symbol | Use case |
|--------|----------|
| `convert_gtf_to_gff3(source, output, ...)` | Convert GTF → GFF3 file |
| `describe_convert_mode(source, ...)` | Mode label + `SourceInfo` for logging |

Converted GFF3 can be passed directly to `compute_gff_stats` for annotrieve-compatible statistics.

## Python library usage

### Quick start

```python
from gffy import compute_gff_stats

stats = compute_gff_stats("/path/to/annotation.gff3.gz")

coding = stats["gene_category_stats"].get("coding", {})
print(coding.get("total_count", 0))

for ttype, tstats in stats["transcript_type_stats"].items():
    print(ttype, tstats["total_count"])
```

### Primary API

```python
compute_gff_stats(
    source: str,
    *,
    force_low_memory: bool = False,
) -> dict
```

- `source` — local file path or `http(s)`/`ftp` URL (plain or gzip).
- `force_low_memory` — always use sort-to-temp + per-seqid processing (see [Memory modes](#memory-modes)).
- **Returns** a dict with `gene_category_stats` and `transcript_type_stats` (annotrieve-compatible).

### Additional exports

| Symbol | Use case |
|--------|----------|
| `compute_gff_stats_from_lines(lines)` | Stats from an in-memory iterable of GFF lines |
| `inspect_source(source)` | Size/metadata for a path or URL (`SourceInfo`) |
| `is_big(info)` | Whether auto low-memory mode applies |
| `describe_compute_mode(source, force_low_memory=...)` | Human-readable mode label + `SourceInfo` |
| `build_sorted_gff(source, output_path, info=None)` | Low-level: sort into a gzipped file you manage |
| `convert_gtf_to_gff3(source, output, ...)` | Convert GTF → GFF3 (see [GTF → GFF3 conversion](#gtf--gff3-conversion)) |
| `describe_convert_mode(source, ...)` | Mode label for conversion |
| `build_custom_annotation(source, ...)` | Annotrieve Favorites `CustomAnnotation` dict (offline import) |
| `derive_custom_name(source, override=None)` | Display name from local path or URL path only (no host) for custom annotations |
| `compute_bulk_stats(sources, output, ...)` | Bulk stats → JSONL (see [Bulk stats from a list](#bulk-stats-from-a-list)) |
| `read_source_list(path_or_dash)` | Parse a newline-separated source list (`SourceListResult`) |
| `warn_source_list_issues(result)` | Print duplicate-line warnings to stderr |
| `ErrorInfo`, `classify_exception` | Structured error classification |
| `__version__` | Package version string |

### Validate output JSON

The output matches [`schema.json`](schema.json) at the repo root. When installed, a copy is shipped inside the package:

```python
from importlib.resources import files
import json

schema = json.loads(
    files("gffy").joinpath("schema.json").read_text(encoding="utf-8")
)
```

### Large files in code

```python
stats = compute_gff_stats(
    "https://example.com/large.gff.gz",
    force_low_memory=True,
)
```

## Memory modes

| Mode | When | Peak RAM |
|------|------|----------|
| **single-pass** | gz ≤ 100 MB and plain ≤ 1 GB (by size check) | Low for small/medium files |
| **low-memory** | Above thresholds, or `force_low_memory=True` / `--low-memory` | Bounded per seqid (~500 MB target) |

Low-memory mode:

1. Sorts the GFF by seqid with system `sort` (`LC_ALL=C`, `-S 200M` buffer cap) into a **temporary** `.gff.gz`.
2. Streams one seqid at a time, updates global stats, frees per-seqid scratch.
3. **Deletes** the temporary sorted file when done (nothing persisted on disk).

Environment overrides:

- `GFFY_GZ_THRESHOLD_BYTES` (default `104857600` — 100 MB)
- `GFFY_PLAIN_THRESHOLD_BYTES` (default `1073741824` — 1 GB)

URL auto-detection uses `HEAD` `Content-Length` when available; missing length defaults to single-pass.

## Output structure

Top-level keys:

- `gene_category_stats` — `coding`, `non_coding`, `pseudogene` (each present only if count &gt; 0)
- `transcript_type_stats` — keyed by transcript type (e.g. `mRNA`), sorted by descending `total_count`

Example (abbreviated):

```json
{
  "gene_category_stats": {
    "coding": {
      "total_count": 22178,
      "length_stats": { "min": 10, "max": 2960899, "mean": 48306.64 },
      "biotype_counts": { "protein_coding": 20000 },
      "transcript_type_counts": { "mRNA": 66153 }
    }
  },
  "transcript_type_stats": {
    "mRNA": {
      "total_count": 66153,
      "length_stats": { "min": 100, "max": 50000, "mean": 3500.0 },
      "biotype_counts": { "protein_coding": 60000 },
      "associated_genes": {
        "total_count": 20000,
        "gene_categories": { "coding": 20000 }
      },
      "exon_stats": { "total_count": 500000, "length": { "min": 10, "max": 5000, "mean": 150.0 } },
      "cds_stats": { "..." : "..." }
    }
  }
}
```

Gene categories:

- **coding** — genes with CDS-bearing transcripts or `protein_coding` biotype
- **non_coding** — genes with exons but no CDS
- **pseudogene** — `pseudogene` feature type

## Development

```bash
pip install -e ".[dev]"
pytest -v
python -m build    # smoke-test packaging
flake8
black .
```

CI runs on push and pull requests to `main` and `master` (Python 3.9 and 3.12).

## Releasing to PyPI

1. Bump `version` in [`pyproject.toml`](pyproject.toml) and [`src/gffy/__init__.py`](src/gffy/__init__.py).
2. Commit, tag (e.g. `v0.1.1`), and push the tag.
3. Create a **GitHub Release** from that tag (publish event triggers [`.github/workflows/publish.yml`](.github/workflows/publish.yml)).

**One-time PyPI trusted publishing setup:**

1. Register the project on [PyPI](https://pypi.org/) as `gffy` (if needed).
2. PyPI → **Publishing** → add a trusted publisher for GitHub:
   - Owner: `guigolab`
   - Repository: `gffy`
   - Workflow: `publish.yml`
   - Environment: `pypi` (matches the workflow `environment.name`)

No long-lived `PYPI_API_TOKEN` is required when trusted publishing is configured.

## License

MIT — see [LICENSE](LICENSE).
