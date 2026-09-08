"""Tests for the combined METEOR_USE_RUST / --use-rust flag and fallback behavior."""

from __future__ import annotations

import json
import lzma
import sys
import tempfile
from pathlib import Path

import pytest

try:
    import meteor_core

    HAS_METEOR_CORE = True
except ImportError:
    meteor_core = None
    HAS_METEOR_CORE = False

import meteor.counter as counter_module
from meteor.counter import Counter
from meteor.session import Component

FIXTURE_DIR = Path(__file__).resolve().parent / "data" / "fixtures"
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


def _run_counter(use_rust_counter: bool) -> dict[int, float]:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
        dir=FIXTURE_DIR, prefix="rust_combined_"
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
        return _parse_tsv(count_file)


def test_use_rust_parser_flags_exist(tmp_path: Path) -> None:
    """The combined --use-rust flag is present on mapping and strain subparsers."""
    from meteor.meteor import get_arguments

    in_dir = tmp_path / "in"
    ref_dir = tmp_path / "ref"
    out_dir = tmp_path / "out"
    for d in (in_dir, ref_dir, out_dir):
        d.mkdir()

    mapping_args = get_arguments(
        [
            "mapping",
            "-i",
            str(in_dir),
            "-r",
            str(ref_dir),
            "-o",
            str(out_dir),
            "--use-rust",
        ]
    )
    assert mapping_args.use_rust is True
    assert mapping_args.use_rust_counter is False

    strain_args = get_arguments(
        [
            "strain",
            "-i",
            str(in_dir),
            "-r",
            str(ref_dir),
            "-o",
            str(out_dir),
            "--use-rust",
        ]
    )
    assert strain_args.use_rust is True
    assert strain_args.use_rust_variant_calling is False


@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="module-level monkeypatching is os-agnostic but keep test simple",
)
def test_counter_fallback_when_meteor_core_missing(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """If meteor_core is unavailable, the Rust counter logs a warning and falls back."""
    real_meteor_core = counter_module.meteor_core
    monkeypatch.setattr(counter_module, "meteor_core", None)
    try:
        python_counts = _run_counter(use_rust_counter=False)
        with caplog.at_level("WARNING", logger="meteor.counter"):
            fallback_counts = _run_counter(use_rust_counter=True)

        assert fallback_counts == python_counts
        assert "not available" in caplog.text
        assert "Falling back to Python" in caplog.text
    finally:
        monkeypatch.setattr(counter_module, "meteor_core", real_meteor_core)


def test_counter_fallback_on_rust_panic(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the Rust counter call raises, the counter falls back to Python."""
    real_meteor_core = counter_module.meteor_core

    class FailingCore:
        @staticmethod
        def count_msp(*_args, **_kwargs) -> None:
            raise RuntimeError("simulated Rust panic")

    monkeypatch.setattr(counter_module, "meteor_core", FailingCore())
    try:
        python_counts = _run_counter(use_rust_counter=False)
        fallback_counts = _run_counter(use_rust_counter=True)
        assert fallback_counts == python_counts
    finally:
        monkeypatch.setattr(counter_module, "meteor_core", real_meteor_core)
