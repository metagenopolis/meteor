"""Tests for lazy CRAM record streaming from Rust."""

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import meteor_core

FIXTURES = Path(__file__).parent / "data" / "fixtures"
CRAM = FIXTURES / "sample.cram"
REF = FIXTURES / "reference.fa"
TIME_BIN = "/usr/bin/time"
EVIDENCE_DIR = (
    Path(__file__).parent.parent
    / ".omo"
    / "evidence"
    / "meteor-rust-acceleration"
    / "streaming"
)


def _streaming_script() -> str:
    return """
import json, sys
sys.path.insert(0, '.')
import meteor_core
it = meteor_core.stream_cram_records(sys.argv[1], sys.argv[2])
count = 0
total_len = 0
total_ops = 0
for rec in it:
    count += 1
    for _, length in rec.cigar_tuples:
        total_len += length
        total_ops += 1
print(json.dumps({"count": count, "total_len": total_len, "total_ops": total_ops}))
"""


def _buffered_script() -> str:
    return """
import json, sys
sys.path.insert(0, '.')
import meteor_core
records = meteor_core.cram_records(sys.argv[1], sys.argv[2])
count = len(records)
total_len = 0
total_ops = 0
for rec in records:
    for _, length in rec.cigar_tuples:
        total_len += length
        total_ops += 1
print(json.dumps({"count": count, "total_len": total_len, "total_ops": total_ops}))
"""


def _run_with_rss(script: str, cram: Path, ref: Path, tmp_path: Path) -> tuple[dict, int]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    script_path = tmp_path / "script.py"
    script_path.write_text(script)
    cmd = [TIME_BIN, "-l", sys.executable, str(script_path), str(cram), str(ref)]
    result = subprocess.run(
        cmd,
        cwd=Path(__file__).parent.parent,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"subprocess failed: {result.stderr}")
    metrics = json.loads(result.stdout.strip().splitlines()[-1])
    match = re.search(r"maximum resident set size\s+(\d+)", result.stderr)
    if not match:
        raise RuntimeError(f"could not parse max RSS from: {result.stderr}")
    max_rss_kb = int(match.group(1)) // 1024  # bytes -> KB
    return metrics, max_rss_kb


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


@pytest.mark.skipif(not shutil.which("/usr/bin/time"), reason="/usr/bin/time not available")
def test_streaming_peak_rss_lower_than_buffered(tmp_path: Path) -> None:
    """Peak RSS of streaming is at most 80% of the buffered path on the same fixture."""
    stream_metrics, stream_rss_kb = _run_with_rss(
        _streaming_script(), CRAM, REF, tmp_path / "stream"
    )
    buffer_metrics, buffer_rss_kb = _run_with_rss(
        _buffered_script(), CRAM, REF, tmp_path / "buffer"
    )

    assert stream_metrics == buffer_metrics

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    evidence = {
        "method": "/usr/bin/time -l (max resident set size in KB)",
        "fixture": str(CRAM),
        "streamed_rss_kb": stream_rss_kb,
        "buffered_rss_kb": buffer_rss_kb,
        "streamed_fraction": stream_rss_kb / buffer_rss_kb if buffer_rss_kb else None,
        "metrics": stream_metrics,
    }
    (EVIDENCE_DIR / "rss_measurement.json").write_text(json.dumps(evidence, indent=2))

    assert stream_rss_kb <= 0.8 * buffer_rss_kb, (
        f"streaming RSS ({stream_rss_kb} KB) is not <= 80% of buffered "
        f"({buffer_rss_kb} KB)"
    )
