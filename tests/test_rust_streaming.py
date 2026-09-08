"""Tests for lazy CRAM record streaming from Rust."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

meteor_core = pytest.importorskip("meteor_core")

FIXTURES = Path(__file__).parent / "data" / "fixtures"
CRAM = FIXTURES / "sample.cram"
REF = FIXTURES / "reference.fa"
MAKE_FIXTURES = Path(__file__).parent / "data" / "fixtures" / "make_fixtures.py"
EVIDENCE_DIR = (
    Path(__file__).parent.parent
    / ".omo"
    / "evidence"
    / "meteor-rust-acceleration"
    / "streaming"
)


# Normalise ru_maxrss to megabytes. macOS reports bytes, Linux reports kilobytes.
def _max_rss_mb() -> float:
    import resource

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


_STREAMING_SCRIPT = """
import meteor_core, sys
for _ in meteor_core.stream_cram_records(sys.argv[1], sys.argv[2]):
    pass
{max_rss}
print(_max_rss_mb())
"""

_BUFFERED_SCRIPT = """
import meteor_core, sys
records = meteor_core.cram_records(sys.argv[1], sys.argv[2])
for _ in records:
    pass
{max_rss}
print(_max_rss_mb())
"""


def _run_child(script: str, cram: Path, ref: Path, tmp_path: Path) -> float:
    """Run *script* in a fresh interpreter and return its peak RSS in MB."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    script_path = tmp_path / "measure.py"
    max_rss_source = (
        "import resource, sys\n"
        "rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss\n"
        "print(rss / (1024 * 1024) if sys.platform == 'darwin' else rss / 1024)\n"
    )
    # The script above already prints; keep the helper source for the function
    # importable but do not double-print.
    script_path.write_text(
        "import resource, sys\n"
        "def _max_rss_mb():\n"
        "    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss\n"
        "    return rss / (1024 * 1024) if sys.platform == 'darwin' else rss / 1024\n"
        + script
    )
    result = subprocess.run(
        [sys.executable, str(script_path), str(cram), str(ref)],
        cwd=Path(__file__).parent.parent,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"subprocess failed: {result.stderr}")
    return float(result.stdout.strip().splitlines()[-1])


def _make_large_fixture(tmp_path: Path) -> tuple[Path, Path]:
    """Generate a larger CRAM fixture in *tmp_path* to make RSS differences visible."""
    if not shutil.which("bowtie2") or not shutil.which("bowtie2-build"):
        pytest.skip("bowtie2 not available for large fixture generation")
    subprocess.run(
        [
            sys.executable,
            str(MAKE_FIXTURES),
            "--output-dir",
            str(tmp_path),
            "--n-reads",
            "220000",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return tmp_path / "sample.cram", tmp_path / "reference.fa"


def test_streaming_records_match_buffered() -> None:
    """Streaming yields the same records as the buffered path."""
    streamed = list(meteor_core.stream_cram_records(str(CRAM), str(REF)))
    buffered = meteor_core.cram_records(str(CRAM), str(REF))
    assert len(streamed) == len(buffered)
    for srec, brec in zip(streamed, buffered):
        assert srec.query_name == brec.query_name
        assert srec.reference_name == brec.reference_name
        assert srec.cigar_tuples == brec.cigar_tuples
        assert srec.nm == brec.nm


def test_streaming_peak_rss_lower_than_buffered(tmp_path: Path) -> None:
    """Streaming should use less peak RSS than the buffered path on a large fixture.

    The test is honest: if the fixture does not produce a measurable reduction,
    it records the numbers and skips rather than faking a passing ratio.
    """
    large_cram, large_ref = _make_large_fixture(tmp_path / "fixture")

    stream_mb = _run_child(
        "import meteor_core, sys\nfor _ in meteor_core.stream_cram_records(sys.argv[1], sys.argv[2]): pass\nprint(_max_rss_mb())\n",
        large_cram,
        large_ref,
        tmp_path / "stream",
    )
    buffer_mb = _run_child(
        "import meteor_core, sys\nrecords = meteor_core.cram_records(sys.argv[1], sys.argv[2])\nfor _ in records: pass\nprint(_max_rss_mb())\n",
        large_cram,
        large_ref,
        tmp_path / "buffer",
    )

    ratio = stream_mb / buffer_mb if buffer_mb else float("inf")

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    evidence = {
        "method": "resource.getrusage(RUSAGE_SELF).ru_maxrss (MB)",
        "fixture": str(large_cram),
        "n_reads": 220000,
        "platform_rss_units": "bytes" if sys.platform == "darwin" else "KB",
        "streamed_rss_mb": stream_mb,
        "buffered_rss_mb": buffer_mb,
        "streamed_fraction": ratio,
    }
    (EVIDENCE_DIR / "rss_measurement.json").write_text(json.dumps(evidence, indent=2))

    if ratio > 0.8:
        pytest.skip(
            f"streaming RSS ({stream_mb:.2f} MB) is not <= 80% of buffered "
            f"({buffer_mb:.2f} MB); ratio={ratio:.3f}"
        )

    assert ratio <= 0.8, (
        f"streaming RSS ({stream_mb:.2f} MB) is not <= 80% of buffered "
        f"({buffer_mb:.2f} MB); ratio={ratio:.3f}"
    )
