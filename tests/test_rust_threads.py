"""Verify METEOR_RUST_THREADS env does not change CRAM decode results."""

from __future__ import annotations

from pathlib import Path

import pytest

meteor_core = pytest.importorskip("meteor_core")

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"


def _count_msp() -> tuple[int, dict[int, float]]:
    """Run the Rust counter on the fixture and return the total + per-gene map."""
    identity_threshold = 0.95
    counting_type = "smart_shared"
    result = meteor_core.count_msp(
        str(CRAM), str(REF), identity_threshold, counting_type
    )
    counts = {record.gene_id: record.count for record in result.gene_counts}
    return result.counted_reads, counts


@pytest.mark.parametrize("thread_env", ["1", "4", "999999", "abc"])
def test_rust_threads_env_parity(monkeypatch: pytest.MonkeyPatch, thread_env: str) -> None:
    """Changing METEOR_RUST_THREADS must not alter decoded counts."""
    # Given: a baseline run with the variable unset.
    monkeypatch.delenv("METEOR_RUST_THREADS", raising=False)
    baseline_reads, baseline_counts = _count_msp()

    # When: the same counter is called under a controlled thread setting.
    monkeypatch.setenv("METEOR_RUST_THREADS", thread_env)
    result_reads, result_counts = _count_msp()

    # Then: the outputs are byte-for-byte identical to the baseline.
    assert result_reads == baseline_reads
    assert result_counts == baseline_counts
