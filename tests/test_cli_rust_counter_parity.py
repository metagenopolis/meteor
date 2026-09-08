"""CLI-level parity test for the METEOR_USE_RUST_COUNTER flag."""

from __future__ import annotations

import json
import lzma
import tempfile
from pathlib import Path

import pytest

from meteor.counter import Counter
from meteor.session import Component

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF_JSON = FIXTURE_DIR / "reference.json"
STAGE1_JSON = FIXTURE_DIR / "sample_census_stage_1.json"


def _parse_tsv(path: Path) -> dict[int, float]:
    counts: dict[int, float] = {}
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            gene_id, _length, value = line.strip().split("\t")
            counts[int(gene_id)] = float(value)
    return counts


def _run_counter(use_rust_counter: bool, identity_threshold: float) -> dict[int, float]:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(dir=FIXTURE_DIR, prefix="rust_counter_cli_") as tmp_raw:
        tmp_dir = Path(tmp_raw)
        count_file = tmp_dir / "sample.tsv.xz"
        stage1_out = tmp_dir / "sample_census_stage_1.json"
        strain_cram = tmp_dir / "strain.cram"

        meteor = Component
        meteor.threads = 1
        meteor.tmp_path = tmp_dir
        meteor.tmp_dir = tmp_dir
        meteor.mapping_dir = FIXTURE_DIR
        meteor.fastq_dir = FIXTURE_DIR
        meteor.ref_dir = FIXTURE_DIR
        meteor.use_rust_counter = use_rust_counter

        counter = Counter(
            meteor,
            "smart_shared",
            "end-to-end",
            80,
            identity_threshold,
            100,
            100,
        )
        counter.identity_threshold = identity_threshold

        counter.launch_counting(
            CRAM,
            strain_cram,
            count_file,
            ref_json,
            stage1_data,
            stage1_out,
        )
        return _parse_tsv(count_file)


def test_rust_counter_flag_parity():
    identity_threshold = 0.95
    python_counts = _run_counter(use_rust_counter=False, identity_threshold=identity_threshold)
    rust_counts = _run_counter(use_rust_counter=True, identity_threshold=identity_threshold)

    assert set(python_counts.keys()) == set(rust_counts.keys())
    for gene_id in python_counts:
        assert python_counts[gene_id] == pytest.approx(
            rust_counts[gene_id], rel=1e-9, abs=1e-9
        )
