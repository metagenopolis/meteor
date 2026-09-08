"""Spike parity test: rust-htslib reads the fixture CRAM the same way pysam does."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

import meteor_core

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"


def _open_pysam():
    return pysam.AlignmentFile(str(CRAM), "rc", reference_filename=str(REF))


def test_version():
    assert meteor_core.__version__ == "0.1.0"


def test_record_count_parity():
    rust_count = meteor_core.count_cram_records(str(CRAM), str(REF))
    with _open_pysam() as handle:
        py_count = sum(1 for _ in handle)
    assert rust_count == py_count


def test_cigar_parity():
    rust_total_len, rust_total_ops = meteor_core.sum_cram_cigar_lengths(
        str(CRAM), str(REF)
    )

    with _open_pysam() as handle:
        reads = list(handle)
    py_total_len = sum(
        length
        for read in reads
        for _op, length in read.cigartuples
    )
    py_total_ops = sum(len(read.cigartuples) for read in reads)

    assert rust_total_len == py_total_len
    assert rust_total_ops == py_total_ops
