"""Parity tests for Rust-accelerated FASTA/faidx and BED chunking helpers."""

from __future__ import annotations

from pathlib import Path

import pytest
from pysam import FastaFile

import meteor_core

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
FASTA = FIXTURE_DIR / "reference.fa"
BED = FIXTURE_DIR / "regions.bed"


def test_load_fasta_matches_pysam():
    rust_seqs = meteor_core.load_fasta(str(FASTA))
    with FastaFile(filename=str(FASTA)) as fasta:
        for name in fasta.references:
            expected = fasta.fetch(name)
            assert rust_seqs[name] == expected


def test_bed_chunks_preserves_rows_and_balances():
    raw_rows = [
        tuple(int(v) for v in line.strip().split("\t")[:3])
        for line in BED.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]

    for num_chunks in (1, 3, 7, 100):
        chunks = meteor_core.bed_chunks(str(BED), num_chunks)
        flattened = [row for chunk in chunks for row in chunk]
        assert flattened == raw_rows
        assert len(chunks) == min(num_chunks, len(raw_rows))
        sizes = [len(chunk) for chunk in chunks]
        assert max(sizes) - min(sizes) <= 1
        assert all(sizes)


def test_bed_chunks_rejects_zero_chunks():
    with pytest.raises(ValueError, match="num_chunks must be > 0"):
        meteor_core.bed_chunks(str(BED), 0)
