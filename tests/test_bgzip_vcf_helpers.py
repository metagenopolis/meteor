"""Tests for Rust VCF/BCF/bgzip output helpers."""

import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

meteor_core = pytest.importorskip("meteor_core")


def _make_records() -> list[meteor_core.VcfRecord]:
    return [
        meteor_core.VcfRecord(
            chrom="1",
            pos=100,
            id="rs1",
            ref_allele="A",
            alt_alleles=["G"],
            qual=30.0,
            filter="PASS",
            info={"RO": "10", "AO": "5"},
        ),
        meteor_core.VcfRecord(
            chrom="1",
            pos=200,
            id=".",
            ref_allele="C",
            alt_alleles=["T", "G"],
            qual=25.5,
            filter="PASS",
            info={},
        ),
        meteor_core.VcfRecord(
            chrom="2",
            pos=50,
            id="rs2",
            ref_allele="AT",
            alt_alleles=["A"],
            qual=-1.0,
            filter=".",
            info={"DP": "20"},
        ),
    ]


def _record_tuple(rec: pysam.VariantRecord) -> tuple:
    return (
        rec.chrom,
        rec.pos,
        rec.id,
        rec.ref,
        tuple(rec.alts) if rec.alts else (),
        rec.qual,
        list(rec.filter.keys()),
        dict(rec.info),
    )


def _records_equal(left: list, right: list) -> bool:
    if len(left) != len(right):
        return False
    for lrec, rrec in zip(left, right):
        if lrec[:7] != rrec[:7]:
            return False
        if dict(lrec[7]) != dict(rrec[7]):
            return False
    return True


def test_write_vcf_text_roundtrip(tmp_path: Path) -> None:
    vcf_path = tmp_path / "out.vcf"
    records = _make_records()
    meteor_core.write_vcf_text(records, str(vcf_path), "sample1")

    with pysam.VariantFile(str(vcf_path)) as vcf:
        read = [_record_tuple(r) for r in vcf]
    expected = [
        ("1", 100, "rs1", "A", ("G",), 30.0, ["PASS"], {"RO": ("10",), "AO": ("5",)}),
        ("1", 200, None, "C", ("T", "G"), 25.5, ["PASS"], {}),
        ("2", 50, "rs2", "AT", ("A",), None, [], {"DP": ("20",)}),
    ]
    assert _records_equal(read, expected)


def test_bgzip_file_and_tabix_index(tmp_path: Path) -> None:
    vcf_path = tmp_path / "out.vcf"
    gz_path = tmp_path / "out.vcf.gz"
    records = _make_records()
    meteor_core.write_vcf_text(records, str(vcf_path), "sample1")
    meteor_core.bgzip_file(str(vcf_path), str(gz_path))
    assert gz_path.exists()

    pysam.tabix_index(str(gz_path), preset="vcf", force=True)
    index_path = tmp_path / "out.vcf.gz.tbi"
    assert index_path.exists()

    with pysam.VariantFile(str(gz_path)) as vcf:
        read = [_record_tuple(r) for r in vcf]
    assert len(read) == len(records)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_write_bcf_roundtrip_and_index(tmp_path: Path) -> None:
    bcf_path = tmp_path / "out.bcf"
    records = _make_records()
    meteor_core.write_bcf(records, str(bcf_path), "sample1")
    assert bcf_path.exists()

    subprocess.run(["bcftools", "index", str(bcf_path)], check=True)
    assert (tmp_path / "out.bcf.csi").exists()

    expected = [
        ("1", 100, "rs1", "A", ("G",), 30.0, ["PASS"], {"RO": ("10",), "AO": ("5",)}),
        ("1", 200, None, "C", ("T", "G"), 25.5, ["PASS"], {}),
        ("2", 50, "rs2", "AT", ("A",), None, [], {"DP": ("20",)}),
    ]
    with pysam.VariantFile(str(bcf_path)) as vcf:
        read = [_record_tuple(r) for r in vcf]
    assert _records_equal(read, expected)


def test_invalid_vcf_record_raises_typed_error() -> None:
    with pytest.raises(meteor_core.VcfError):
        meteor_core.VcfRecord(
            chrom="1",
            pos=10,
            id=".",
            ref_allele="",
            alt_alleles=["G"],
            qual=10.0,
            filter="PASS",
            info={},
        )
