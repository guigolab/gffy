"""
Sort helpers for large GFF3 files.

Big files are sorted by seqid (system ``sort``) into a temporary gzipped file
for per-seqid streaming; callers remove the temp dir when done.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO
from urllib.request import Request, urlopen

from ._io import is_gzipped_path, is_url


def urlopen_request(
    url: str,
    *,
    method: str = "GET",
    timeout: int = 300,
) -> BinaryIO:
    """Open a URL with User-Agent and timeout; re-raises on failure."""
    from . import __version__

    req = Request(url, method=method)
    req.add_header("User-Agent", f"gffy/{__version__}")
    return urlopen(req, timeout=timeout)  # type: ignore[return-value]

GFFY_GZ_THRESHOLD_BYTES = int(
    os.environ.get("GFFY_GZ_THRESHOLD_BYTES", str(100 * 1024 * 1024))
)
GFFY_PLAIN_THRESHOLD_BYTES = int(
    os.environ.get("GFFY_PLAIN_THRESHOLD_BYTES", str(1024 * 1024 * 1024))
)


@dataclass(frozen=True)
class SourceInfo:
    source: str
    is_url: bool
    size: int | None
    is_gzipped: bool
    etag: str | None = None
    last_modified: str | None = None
    local_path: str | None = None
    mtime_ns: int | None = None


def _head_url(source: str) -> tuple[int | None, str | None, str | None]:
    try:
        with urlopen_request(source, method="HEAD", timeout=60) as resp:
            length = resp.headers.get("Content-Length")
            size = int(length) if length and length.isdigit() else None
            return size, resp.headers.get("ETag"), resp.headers.get("Last-Modified")
    except Exception:
        return None, None, None


def inspect_source(source: str) -> SourceInfo:
    """Gather size and identity metadata for big-file detection."""
    is_gz = is_gzipped_path(source)

    if is_url(source):
        size, etag, last_modified = _head_url(source)
        return SourceInfo(
            source=source,
            is_url=True,
            size=size,
            is_gzipped=is_gz,
            etag=etag,
            last_modified=last_modified,
        )

    path = Path(source).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"GFF file not found: {source}")
    st = path.stat()
    return SourceInfo(
        source=str(path),
        is_url=False,
        size=st.st_size,
        is_gzipped=is_gz or is_gzipped_path(path),
        local_path=str(path),
        mtime_ns=st.st_mtime_ns,
    )


def is_big(info: SourceInfo) -> bool:
    """Return True when the source exceeds configured size thresholds."""
    if info.size is None:
        return False
    if info.is_gzipped:
        return info.size > GFFY_GZ_THRESHOLD_BYTES
    return info.size > GFFY_PLAIN_THRESHOLD_BYTES


def _download_url(source: str, dest: Path) -> None:
    with urlopen_request(source, timeout=300) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def _sort_to_gzip(input_path: Path, output_path: Path, is_gzipped: bool) -> None:
    sort_bin = shutil.which("sort")
    if not sort_bin:
        raise RuntimeError(
            "system 'sort' is required for low-memory mode on large files; "
            "install coreutils or use a smaller file"
        )

    nproc = os.cpu_count() or 1
    tmp_dir = output_path.parent / "sort_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    sort_part = (
        f"LC_ALL=C {shlex.quote(sort_bin)} -t$'\\t' -k1,1 -s -S 200M "
        f"--parallel={nproc} -T {shlex.quote(str(tmp_dir))}"
    )

    if is_gzipped:
        decompress = f"gzip -dc {shlex.quote(str(input_path))}"
        cmd = f"{decompress} | {sort_part} | gzip > {shlex.quote(str(output_path))}"
    else:
        cmd = (
            f"{sort_part} {shlex.quote(str(input_path))} "
            f"| gzip > {shlex.quote(str(output_path))}"
        )

    proc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or "sort pipeline failed")

    try:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    except OSError:
        pass


def build_sorted_gff(source: str, output_path: Path, info: SourceInfo | None = None) -> Path:
    """
    Sort *source* by seqid and write gzipped output to *output_path*.

    *output_path* must live in a temporary directory removed by the caller after use.
    """
    info = info or inspect_source(source)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if info.is_url:
        raw = output_path.parent / (
            f"raw{'.gff.gz' if info.is_gzipped else '.gff'}"
        )
        _download_url(source, raw)
        _sort_to_gzip(raw, output_path, info.is_gzipped)
        try:
            raw.unlink()
        except OSError:
            pass
    else:
        inp = Path(info.local_path or source)
        _sort_to_gzip(inp, output_path, info.is_gzipped or is_gzipped_path(inp))
    return output_path
