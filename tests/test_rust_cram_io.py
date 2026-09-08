"""Parity test for Rust CRAM I/O against pysam."""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

meteor_core = pytest.importorskip("meteor_core")

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"


def test_record_count_and_cigar_parity():
    rust_records = meteor_core.cram_records(str(CRAM), str(REF))

    with pysam.AlignmentFile(str(CRAM), "rc", reference_filename=str(REF)) as handle:
        py_records = list(handle)

    assert len(rust_records) == len(py_records)

    rust_total_len = sum(
        length for record in rust_records for _op, length in record.cigar_tuples
    )
    rust_total_ops = sum(len(record.cigar_tuples) for record in rust_records)
    py_total_len = sum(
        length for read in py_records for _op, length in read.cigartuples
    )
    py_total_ops = sum(len(read.cigartuples) for read in py_records)

    assert rust_total_len == py_total_len
    assert rust_total_ops == py_total_ops


def test_per_record_fields_parity():
    rust_records = meteor_core.cram_records(str(CRAM), str(REF))

    with pysam.AlignmentFile(str(CRAM), "rc", reference_filename=str(REF)) as handle:
        py_records = list(handle)

    assert len(rust_records) == len(py_records)

    for rust_record, py_record in zip(rust_records, py_records):
        assert rust_record.query_name == py_record.query_name
        assert rust_record.reference_name == py_record.reference_name
        assert rust_record.cigar_tuples == list(py_record.cigartuples)
        rust_nm = rust_record.nm
        py_nm = py_record.get_tag("NM") if py_record.has_tag("NM") else None
        assert rust_nm == py_nm
