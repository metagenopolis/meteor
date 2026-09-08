"""Parity tests for the Rust region-parallel freebayes dispatcher."""

import shutil
import subprocess
from pathlib import Path

import pysam
import pytest

meteor_core = pytest.importorskip("meteor_core")

FIXTURES = Path(__file__).parent / "data" / "fixtures"
CRAM = FIXTURES / "sample.cram"
REF = FIXTURES / "reference.fa"
BED = FIXTURES / "regions.bed"
FREEBAYES = "freebayes"
BCFTOOLS = "bcftools"
TABIX = "tabix"


def _index_reference(ref_path: Path) -> None:
    """Build a FASTA index if it is missing."""
    fai = Path(f"{ref_path}.fai")
    if not fai.exists():
        pysam.faidx(str(ref_path))


def _run_freebayes_serial(
    cram_path: Path,
    ref_path: Path,
    bed_path: Path,
    output_path: Path,
    min_snp_depth: int = 1,
    min_frequency: float = 0.1,
    ploidy: int = 1,
) -> None:
    """Run a single freebayes process over the whole BED for reference output."""
    cmd = [
        FREEBAYES,
        "-i",
        "-u",
        "--pooled-continuous",
        "--haplotype-length",
        "0",
        "--min-alternate-count",
        "1",
        "--min-coverage",
        str(min_snp_depth),
        "--min-alternate-fraction",
        str(min_frequency),
        "--min-mapping-quality",
        "0",
        "--use-duplicate-reads",
        "-t",
        str(bed_path),
        "-p",
        str(ploidy),
        "-f",
        str(ref_path),
        "-b",
        str(cram_path),
    ]
    with output_path.open("wb") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"freebayes serial failed (exit {result.returncode}): "
            f"{result.stderr.decode('utf-8', errors='replace')}"
        )


def _normalize_vcf(input_vcf: Path, ref_path: Path, output_vcf: Path) -> None:
    """Normalize and sort a VCF so records can be compared."""
    norm_cmd = [
        BCFTOOLS,
        "norm",
        "-f",
        str(ref_path),
        str(input_vcf),
    ]
    sort_cmd = [
        BCFTOOLS,
        "sort",
        "-Oz",
        "-o",
        str(output_vcf),
        "-T",
        str(output_vcf.parent),
        "-",
    ]
    norm = subprocess.run(norm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if norm.returncode != 0:
        raise RuntimeError(
            f"bcftools norm failed (exit {norm.returncode}): "
            f"{norm.stderr.decode('utf-8', errors='replace')}"
        )
    sort = subprocess.run(sort_cmd, input=norm.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if sort.returncode != 0:
        raise RuntimeError(
            f"bcftools sort failed (exit {sort.returncode}): "
            f"{sort.stderr.decode('utf-8', errors='replace')}"
        )
    subprocess.run([TABIX, "-p", "vcf", str(output_vcf)], check=True)


def _vcf_records(path: Path):
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            yield (
                rec.chrom,
                rec.pos,
                rec.ref,
                tuple(rec.alts) if rec.alts else (),
                rec.qual,
            )


@pytest.mark.skipif(
    not shutil.which(FREEBAYES) or not shutil.which(BCFTOOLS),
    reason="freebayes/bcftools not available",
)
def test_call_variants_parallel_matches_serial(tmp_path: Path) -> None:
    """Rust parallel dispatcher output matches a single freebayes invocation after norm+sort."""
    import shutil

    _index_reference(REF)

    serial_vcf = tmp_path / "serial.vcf"
    rust_vcf = tmp_path / "rust.vcf"

    _run_freebayes_serial(CRAM, REF, BED, serial_vcf)

    options = meteor_core.FreebayesOptions(
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
    )
    meteor_core.call_variants_parallel(
        str(CRAM),
        str(REF),
        str(BED),
        FREEBAYES,
        options,
        n_threads=2,
        output_path=str(rust_vcf),
    )

    serial_norm = tmp_path / "serial.norm.vcf.gz"
    rust_norm = tmp_path / "rust.norm.vcf.gz"
    _normalize_vcf(serial_vcf, REF, serial_norm)
    _normalize_vcf(rust_vcf, REF, rust_norm)

    serial_records = list(_vcf_records(serial_norm))
    rust_records = list(_vcf_records(rust_norm))
    assert len(serial_records) == len(rust_records)
    assert serial_records == rust_records


def test_call_variants_parallel_reports_failure() -> None:
    """A failing freebayes chunk is reported as a typed error without hanging."""
    options = meteor_core.FreebayesOptions(
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
    )
    with pytest.raises(meteor_core.FreebayesError):
        meteor_core.call_variants_parallel(
            str(CRAM),
            str(REF),
            str(BED),
            "/nonexistent/freebayes/binary",
            options,
            n_threads=2,
            output_path=None,
        )
