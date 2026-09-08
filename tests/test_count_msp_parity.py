"""Parity test for the Rust MSP/gene counting hot loop against meteor Counter."""

from __future__ import annotations

import json
import lzma
import tempfile
from pathlib import Path

import pytest

meteor_core = pytest.importorskip("meteor_core")
from meteor.counter import Counter
from meteor.session import Component

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"
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


def _run_python_counter(
    identity_threshold: float, counting_type: str
) -> dict[int, float]:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="bench_counter_test_"
    ) as tmp_raw:
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

        counter = Counter(
            meteor,
            counting_type,
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


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_parity(counting_type: str) -> None:
    identity_threshold = 0.95
    expected = _run_python_counter(identity_threshold, counting_type)
    result = meteor_core.count_msp(
        str(CRAM), str(REF), identity_threshold, counting_type
    )

    expected_by_gene = {gene_id: value for gene_id, value in expected.items()}
    for record in result.gene_counts:
        assert record.gene_id in expected_by_gene
        assert record.count == pytest.approx(
            expected_by_gene[record.gene_id], rel=1e-9, abs=1e-9
        )
    assert len(result.gene_counts) == len(expected_by_gene)
