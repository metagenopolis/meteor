"""End-to-end parity regression suite for Rust-accelerated paths.

Runs the same Meteor operation through the Python and Rust implementations and
asserts that the final outputs are equivalent.
"""

from __future__ import annotations

import json
import lzma
import shutil
import tempfile
from pathlib import Path

import pytest

meteor_core = pytest.importorskip("meteor_core")
from meteor.counter import Counter
from meteor.session import Component
from meteor.variantcalling import VariantCalling

FIXTURES = Path(__file__).parent / "data" / "fixtures"
CRAM = FIXTURES / "sample.cram"
REF_JSON = FIXTURES / "reference.json"
STAGE1_JSON = FIXTURES / "sample_census_stage_1.json"


def _parse_counts(path: Path) -> dict[int, float]:
    counts: dict[int, float] = {}
    opener = lzma.open if path.suffix == ".xz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            gene_id, _length, value = line.strip().split("\t")
            counts[int(gene_id)] = float(value)
    return counts


def _run_counter(use_rust_counter: bool) -> tuple[dict[int, float], dict]:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))
    tmp_dir = Path(tempfile.mkdtemp(dir=FIXTURES, prefix="parity_counter_"))
    try:
        count_file = tmp_dir / "sample.tsv.xz"
        stage1_out = tmp_dir / "sample_census_stage_1.json"
        strain_cram = tmp_dir / "strain.cram"

        meteor = Component(
            threads=1,
            tmp_path=tmp_dir,
            tmp_dir=tmp_dir,
            mapping_dir=FIXTURES,
            fastq_dir=FIXTURES,
            ref_dir=FIXTURES,
        )
        meteor.use_rust_counter = use_rust_counter

        counter = Counter(
            meteor,
            "smart_shared",
            "end-to-end",
            80,
            None,
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
        return (
            _parse_counts(count_file),
            json.loads(stage1_out.read_text(encoding="utf-8")),
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _run_strain(use_rust_variant_calling: bool) -> tuple[Path, Path]:
    ref_json = json.loads(REF_JSON.read_text(encoding="utf-8"))
    stage1_data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))
    sample_name = stage1_data["sample_info"]["sample_name"]
    tmp_dir = Path(tempfile.mkdtemp(dir=FIXTURES, prefix="parity_strain_"))
    out_dir = tmp_dir / "out"
    out_dir.mkdir()

    meteor = Component(
        threads=1,
        tmp_path=tmp_dir,
        tmp_dir=tmp_dir,
        ref_dir=FIXTURES,
        strain_dir=out_dir,
    )
    meteor.use_rust_variant_calling = use_rust_variant_calling

    census = {
        "mapped_sample_dir": FIXTURES,
        "directory": out_dir,
        "Stage3FileName": out_dir / f"{sample_name}_census_stage_3.json",
        "census": stage1_data,
        "reference": ref_json,
    }
    variant_caller = VariantCalling(
        meteor,
        census,
        max_depth=100,
        min_depth=3,
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
        core_size=100,
    )
    variant_caller.execute()
    return tmp_dir, out_dir / f"{sample_name}.vcf.gz", out_dir / f"{sample_name}_consensus.fasta.xz"


def _read_consensus(path: Path) -> str:
    with lzma.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def _compare_vcfs(left: Path, right: Path) -> None:
    """Compare two VCFs by normalized records using pysam."""
    import pysam

    def _records(path: Path):
        with pysam.VariantFile(str(path)) as vcf:
            for rec in vcf:
                yield (
                    rec.chrom,
                    rec.pos,
                    rec.ref,
                    tuple(rec.alts) if rec.alts else (),
                    round(rec.qual, 3) if rec.qual is not None else None,
                )

    left_records = list(_records(left))
    right_records = list(_records(right))
    assert len(left_records) == len(right_records), (
        f"VCF record counts differ: {len(left_records)} vs {len(right_records)}"
    )
    for l_rec, r_rec in zip(left_records, right_records):
        assert l_rec == r_rec, f"VCF records differ: {l_rec} vs {r_rec}"


def test_counter_python_rust_parity() -> None:
    """Counter output must be identical between Python and Rust paths."""
    python_counts, python_stage1 = _run_counter(use_rust_counter=False)
    rust_counts, rust_stage1 = _run_counter(use_rust_counter=True)

    assert python_counts == pytest.approx(rust_counts, rel=1e-9, abs=1e-9)
    assert python_stage1 == rust_stage1


@pytest.mark.skipif(
    not shutil.which("freebayes") or not shutil.which("bcftools"),
    reason="freebayes/bcftools required for strain parity",
)
def test_strain_python_rust_parity() -> None:
    """Variant-calling outputs must be identical between Python and Rust paths."""
    python_tmp, python_vcf, python_consensus = _run_strain(
        use_rust_variant_calling=False
    )
    rust_tmp, rust_vcf, rust_consensus = _run_strain(use_rust_variant_calling=True)
    try:
        _compare_vcfs(python_vcf, rust_vcf)
        assert _read_consensus(python_consensus) == _read_consensus(rust_consensus)
    finally:
        shutil.rmtree(python_tmp, ignore_errors=True)
        shutil.rmtree(rust_tmp, ignore_errors=True)
