"""Real-data parity tests for phase-2 Rust acceleration.

These tests are gated by the METEOR_REAL_BENCH_DATA environment variable. When
unset, the module is skipped so CI stays green. When set, the placeholder
tests below document the exact parity assertions that will be implemented in
Todo 8.
"""

from __future__ import annotations

import os

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("METEOR_REAL_BENCH_DATA") != "1",
    reason="real-data benchmarks not enabled",
)


def test_counter_tsv_byte_identity() -> None:
    """Placeholder: assert Python and Rust counter TSV outputs are byte-identical."""
    pytest.skip("Will assert counter TSV byte identity between Python and Rust paths.")


def test_vcf_record_equality_after_normalization() -> None:
    """Placeholder: assert VCF records match after bcftools norm -f ref | bcftools sort."""
    pytest.skip(
        "Will assert VCF record equality after bcftools norm -f reference.fa | bcftools sort."
    )


def test_consensus_fasta_identity() -> None:
    """Placeholder: assert consensus FASTA sequences are identical."""
    pytest.skip("Will assert consensus FASTA identity between Python and Rust paths.")


def test_depth_array_equality() -> None:
    """Placeholder: assert per-gene depth arrays are identical."""
    pytest.skip(
        "Will assert per-gene depth array equality between Python and Rust paths."
    )
