"""End-to-end parity test for Rust-accelerated variant calling helpers."""

from __future__ import annotations

import json
import lzma
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

meteor_core = pytest.importorskip("meteor_core")
from meteor.session import Component
from meteor.variantcalling import VariantCalling

DATA_DIR = (
    Path(__file__).resolve().parent.parent / "meteor" / "tests" / "test_variantcalling"
)
EVA71_DIR = DATA_DIR / "eva71"
EVA71_BENCH = DATA_DIR / "eva71_bench"


def _build_vc(tmp_path: Path, use_rust_variant_calling: bool) -> VariantCalling:
    meteor = Component
    meteor.ref_dir = EVA71_DIR
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_path
    meteor.use_rust_variant_calling = use_rust_variant_calling
    meteor.DEFAULT_GAP_CHAR = "?"

    ref_json = json.loads(
        (EVA71_DIR / "eva71_reference.json").read_text(encoding="utf-8")
    )
    census_json = json.loads(
        (EVA71_BENCH / "eva71_bench_census_stage_1.json").read_text(encoding="utf-8")
    )
    sample_info = census_json["sample_info"]
    stage3_dir = tmp_path / sample_info["sample_name"]
    stage3_dir.mkdir(exist_ok=True, parents=True)
    data_dict = {
        "mapped_sample_dir": EVA71_BENCH,
        "census": census_json,
        "directory": stage3_dir,
        "Stage3FileName": stage3_dir / "eva71_bench_census_stage_1.json",
        "reference": ref_json,
    }
    vc = VariantCalling(meteor, data_dict, 100, 3, 1, 0.01, 1, 100)
    vc.matrix_file = EVA71_BENCH / "eva71_bench.tsv.xz"
    return vc


def test_filter_low_cov_sites_parity() -> None:
    cram_file = EVA71_BENCH / "eva71_bench.cram"
    reference_file = EVA71_DIR / "fasta" / "eva71.fasta.gz"

    with tempfile.TemporaryDirectory() as tmp_raw:
        tmp = Path(tmp_raw)
        vc_py = _build_vc(tmp / "py", use_rust_variant_calling=False)
        vc_rust = _build_vc(tmp / "rust", use_rust_variant_calling=True)

        py_low_cov, _ = vc_py.filter_low_cov_sites(cram_file, reference_file)
        rust_low_cov, _ = vc_rust.filter_low_cov_sites(cram_file, reference_file)

        assert_frame_equal(
            py_low_cov.sort_values(by=list(py_low_cov.columns)).reset_index(drop=True),
            rust_low_cov.sort_values(by=list(rust_low_cov.columns)).reset_index(
                drop=True
            ),
            check_like=True,
        )


def test_create_consensus_parity() -> None:
    reference_file = EVA71_DIR / "fasta" / "eva71.fasta.gz"
    vcf_file = EVA71_BENCH / "eva71_bench.vcf.gz"
    bed_file = EVA71_DIR / "database" / "eva71.bed"

    with tempfile.TemporaryDirectory() as tmp_raw:
        tmp = Path(tmp_raw)
        vc_py = _build_vc(tmp / "py", use_rust_variant_calling=False)
        vc_rust = _build_vc(tmp / "rust", use_rust_variant_calling=True)

        low_cov_sites = pd.read_table(
            DATA_DIR / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
        ).set_index("gene_id")

        py_consensus = tmp / "py_consensus.fasta.xz"
        rust_consensus = tmp / "rust_consensus.fasta.xz"

        vc_py.create_consensus(
            reference_file,
            py_consensus,
            low_cov_sites,
            pd.DataFrame(),
            vcf_file,
            bed_file,
        )
        vc_rust.create_consensus(
            reference_file,
            rust_consensus,
            low_cov_sites,
            pd.DataFrame(),
            vcf_file,
            bed_file,
        )

        with lzma.open(py_consensus, "rt", encoding="utf-8") as fh:
            py_text = fh.read()
        with lzma.open(rust_consensus, "rt", encoding="utf-8") as fh:
            rust_text = fh.read()
        assert py_text == rust_text
