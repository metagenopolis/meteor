"""Parity and byte-identity tests for the aggregates-only Rust counter API."""

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
MSP_MAP = FIXTURE_DIR / "msp_map.tsv"


def _msp_map() -> dict[int, str]:
    mapping: dict[int, str] = {}
    with Path.open(MSP_MAP, "rt", encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            msp_name, gene_id, _gene_name, _category = line.strip().split("\t")
            mapping[int(gene_id)] = msp_name
    return mapping


def _normalise(value: float) -> int | float:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def _run_python_counter_bytes(
    identity_threshold: float, counting_type: str
) -> tuple[bytes, dict[int, int]]:
    """Run the Python counter path and return the TSV bytes plus gene lengths."""
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
        lengths: dict[int, int] = {}
        with lzma.open(count_file, "rt", encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                gene_id, gene_length, _value = line.strip().split("\t")
                lengths[int(gene_id)] = int(gene_length)
        return lzma.open(count_file, "rb").read(), lengths


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_aggregates_parity(counting_type: str) -> None:
    identity_threshold = 0.95
    aggregates = meteor_core.count_msp_aggregates(
        str(CRAM), str(MSP_MAP), identity_threshold, counting_type
    )
    reference = meteor_core.count_msp(
        str(CRAM), str(REF), identity_threshold, counting_type
    )

    msp_map = _msp_map()
    ref_by_gene = {
        gc.gene_id: (gc.gene_length, gc.count) for gc in reference.gene_counts
    }

    assert len(aggregates) == len(ref_by_gene)
    for row in aggregates:
        gene_id = int(row.gene)
        assert gene_id in ref_by_gene
        assert row.msp == msp_map[gene_id]
        assert row.count == pytest.approx(
            ref_by_gene[gene_id][1], rel=1e-9, abs=1e-9
        )


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_aggregates_tsv_byte_identical(counting_type: str) -> None:
    identity_threshold = 0.95
    python_tsv, lengths = _run_python_counter_bytes(
        identity_threshold, counting_type
    )
    aggregates = meteor_core.count_msp_aggregates(
        str(CRAM), str(MSP_MAP), identity_threshold, counting_type
    )

    lines = ["gene_id\tgene_length\tvalue\n"]
    for row in sorted(aggregates, key=lambda r: int(r.gene)):
        gene_id = int(row.gene)
        value = _normalise(row.count)
        lines.append(f"{gene_id}\t{lengths[gene_id]}\t{value}\n")
    rust_tsv = "".join(lines).encode("utf-8")

    assert rust_tsv == python_tsv


def _run_counter_with_flag(use_rust_counter: bool) -> bytes:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="rust_agg_fallback_"
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
        meteor.use_rust_counter = use_rust_counter

        counter = Counter(
            meteor,
            "smart_shared",
            "end-to-end",
            80,
            0.95,
            100,
            100,
        )
        counter.identity_threshold = 0.95

        counter.launch_counting(
            CRAM,
            strain_cram,
            count_file,
            ref_json,
            stage1_data,
            stage1_out,
        )
        return lzma.open(count_file, "rb").read()


def test_count_msp_aggregates_stale_extension_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the installed meteor_core lacks count_msp_aggregates, fall back to count_msp."""
    monkeypatch.delattr(meteor_core, "count_msp_aggregates", raising=False)
    python_tsv = _run_counter_with_flag(use_rust_counter=False)
    rust_tsv = _run_counter_with_flag(use_rust_counter=True)
    assert rust_tsv == python_tsv


def test_count_msp_aggregates_bad_counting_type() -> None:
    with pytest.raises(ValueError, match="not a valid counting type"):
        meteor_core.count_msp_aggregates(
            str(CRAM), str(MSP_MAP), 0.95, "unknown_mode"
        )


def test_count_msp_aggregates_missing_cram() -> None:
    with pytest.raises(OSError):
        meteor_core.count_msp_aggregates(
            str(FIXTURE_DIR / "does_not_exist.cram"),
            str(MSP_MAP),
            0.95,
            "smart_shared",
        )
