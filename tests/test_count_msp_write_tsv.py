"""Parity tests for the Rust TSV-writing counter API."""

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


def _run_python_counter_bytes(
    identity_threshold: float, counting_type: str
) -> tuple[bytes, dict[int, int]]:
    """Run the Python counter path and return the raw TSV bytes plus gene lengths."""
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
        return count_file.read_bytes(), lengths


def _run_rust_write_tsv(
    identity_threshold: float, counting_type: str
) -> tuple[bytes, int, int]:
    """Run the Rust TSV writer and return raw bytes, row count and counted reads."""
    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="rust_tsv_test_"
    ) as tmp_raw:
        tmp_dir = Path(tmp_raw)
        out_tsv = tmp_dir / "sample.tsv.xz"
        row_count, counted_reads = meteor_core.count_msp_write_tsv(
            str(CRAM), str(MSP_MAP), str(out_tsv), identity_threshold, counting_type
        )
        return out_tsv.read_bytes(), row_count, counted_reads


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_write_tsv_byte_identical(counting_type: str) -> None:
    """Rust-written TSV must match the Python Counter.write_stat output byte-for-byte."""
    identity_threshold = 0.95
    python_tsv, _lengths = _run_python_counter_bytes(identity_threshold, counting_type)
    rust_tsv, _row_count, _counted_reads = _run_rust_write_tsv(
        identity_threshold, counting_type
    )
    assert rust_tsv == python_tsv


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_write_tsv_row_count(counting_type: str) -> None:
    """The returned row count equals the number of gene data rows written."""
    identity_threshold = 0.95
    rust_tsv, row_count, _counted_reads = _run_rust_write_tsv(
        identity_threshold, counting_type
    )
    lines = lzma.decompress(rust_tsv).decode("utf-8").splitlines()
    assert lines[0] == "gene_id\tgene_length\tvalue"
    assert row_count == len(lines) - 1


@pytest.mark.parametrize("counting_type", ["smart_shared", "unique", "total"])
def test_count_msp_write_tsv_counted_reads_match_python(counting_type: str) -> None:
    """The Rust counted_reads must equal the Python path's counted_reads."""
    identity_threshold = 0.95
    _python_tsv, _lengths = _run_python_counter_bytes(identity_threshold, counting_type)
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="rust_tsv_reads_parity_"
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
        meteor.use_rust_counter = True

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
        rust_census = json.loads(stage1_out.read_text(encoding="utf-8"))

    python_tmp = Path(tempfile.mkdtemp(dir=FIXTURE_DIR, prefix="python_reads_parity_"))
    try:
        python_count_file = python_tmp / "sample.tsv.xz"
        python_stage1_out = python_tmp / "sample_census_stage_1.json"
        python_strain_cram = python_tmp / "strain.cram"
        meteor = Component
        meteor.threads = 1
        meteor.tmp_path = python_tmp
        meteor.tmp_dir = python_tmp
        meteor.mapping_dir = FIXTURE_DIR
        meteor.fastq_dir = FIXTURE_DIR
        meteor.ref_dir = FIXTURE_DIR
        meteor.use_rust_counter = False

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
            python_strain_cram,
            python_count_file,
            ref_json,
            stage1_data,
            python_stage1_out,
        )
        python_census = json.loads(python_stage1_out.read_text(encoding="utf-8"))
    finally:
        import shutil

        shutil.rmtree(python_tmp, ignore_errors=True)

    assert (
        rust_census["counting"]["counted_reads"]
        == python_census["counting"]["counted_reads"]
    )


def _run_counter_with_flag(use_rust_counter: bool) -> bytes:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="rust_tsv_fallback_"
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


def test_count_msp_write_tsv_stale_extension_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the installed meteor_core lacks count_msp_write_tsv, fall back to aggregates."""
    monkeypatch.delattr(meteor_core, "count_msp_write_tsv", raising=False)
    python_tsv = _run_counter_with_flag(use_rust_counter=False)
    rust_tsv = _run_counter_with_flag(use_rust_counter=True)
    assert rust_tsv == python_tsv


def test_count_msp_write_tsv_unwritable_output(tmp_path: Path) -> None:
    """An unwritable output path must raise a clear OSError from Rust."""
    out_tsv = tmp_path / "does_not_exist" / "out.tsv.xz"
    with pytest.raises(OSError):
        meteor_core.count_msp_write_tsv(
            str(CRAM), str(MSP_MAP), str(out_tsv), 0.95, "smart_shared"
        )


def test_count_msp_write_tsv_bad_counting_type() -> None:
    """A bad counting type must raise a clear ValueError."""
    with pytest.raises(ValueError, match="not a valid counting type"):
        meteor_core.count_msp_write_tsv(
            str(CRAM), str(MSP_MAP), "/dev/null", 0.95, "unknown_mode"
        )
