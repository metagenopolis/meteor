"""Real-data parity tests for phase-2 Rust acceleration.

These tests are gated by the METEOR_REAL_BENCH_DATA environment variable. When
unset, the module is skipped so CI stays green. When set, the suite compares the
Python and Rust user-facing paths on the real benchmark CRAM/reference data.
"""

from __future__ import annotations

import csv
import difflib
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
import pysam
import pytest

from meteor.counter import Counter
from meteor.session import Component
from meteor.variantcalling import VariantCalling

try:
    import meteor_core
except ImportError:
    meteor_core = None  # type: ignore[assignment]

pytestmark = pytest.mark.skipif(
    os.environ.get("METEOR_REAL_BENCH_DATA") != "1",
    reason="real-data benchmarks not enabled",
)

REQUIRED_ENVS = (
    "METEOR_BENCH_CRAM",
    "METEOR_BENCH_REF",
    "METEOR_BENCH_CATALOGUE",
    "METEOR_BENCH_MSP_MAP",
)

REPO_ROOT = Path(__file__).resolve().parent.parent
EVIDENCE_DIR = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-speedup-phase2" / "realdata"
)

REPORT_DIR = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-speedup-phase2" / "task-8-parity"
)

FREEBAYES = "freebayes"
BCFTOOLS = "bcftools"
TABIX = "tabix"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _env_paths() -> dict[str, Path]:
    """Return the required benchmark paths, skipping if any are missing."""
    missing = [name for name in REQUIRED_ENVS if not os.environ.get(name)]
    if missing:
        pytest.skip(f"missing environment variable(s): {', '.join(missing)}")
    return {name: Path(os.environ[name]) for name in REQUIRED_ENVS}


def _evidence_dir() -> Path:
    """Ensure the evidence directory exists and return it."""
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    return EVIDENCE_DIR


def _write_text_diff(
    path: Path, left_label: str, left: str, right_label: str, right: str
) -> None:
    """Write a unified text diff between *left* and *right* to *path*."""
    diff = difflib.unified_diff(
        left.splitlines(keepends=True),
        right.splitlines(keepends=True),
        fromfile=left_label,
        tofile=right_label,
    )
    path.write_text("".join(diff), encoding="utf-8")


def _stage1_and_ref(
    mapping_dir: Path, ref_dir: Path
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    """Locate the stage-1 JSON and reference JSON for the benchmark sample."""
    stage1_json = next(mapping_dir.rglob("*_census_stage_1.json"))
    stage1_data = json.loads(stage1_json.read_text(encoding="utf-8"))
    ref_json_path = next(ref_dir.glob("*_reference.json"))
    ref_json = json.loads(ref_json_path.read_text(encoding="utf-8"))
    return stage1_json, stage1_data, ref_json


def _reference_fasta(ref_dir: Path, tmp_path: Path) -> Path:
    """Return an uncompressed, indexed reference FASTA to use for the tests.

    The reference is never copied into the repository evidence directory; if a
    compressed source has to be decompressed, it is written to *tmp_path*.
    """
    import gzip

    # Prefer a literal 'reference.fa' entry (used by the ref_for_rust layout).
    candidate = ref_dir / "reference.fa"
    if candidate.is_symlink() or candidate.is_file():
        src = candidate
    else:
        src = next(ref_dir.rglob("*.fa*"))

    src = src.resolve()
    out_fa = tmp_path / "reference.fa"
    out_fai = Path(f"{out_fa}.fai")

    with src.open("rb") as fh:
        magic = fh.read(2)
    is_compressed = (
        src.suffix == ".gz"
        or str(src).endswith(".fasta.gz")
        or magic == b"\x1f\x8b"
    )
    if is_compressed:
        with gzip.open(src, "rb") as reader:
            out_fa.write_bytes(reader.read())
    else:
        out_fa.write_bytes(src.read_bytes())

    if not out_fai.is_file():
        pysam.faidx(str(out_fa))
    return out_fa


def _read_bed_intervals(
    catalogue: Path, limit: int | None = None
) -> list[tuple[str, int, int]]:
    """Read (gene_id, start, end) intervals from the catalogue BED."""
    intervals: list[tuple[str, int, int]] = []
    with catalogue.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3:
                continue
            intervals.append((row[0], int(row[1]), int(row[2])))
    if limit is not None:
        intervals = intervals[:limit]
    return intervals


def _write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    """Write intervals as a tab-separated BED file."""
    with path.open("w") as fh:
        for gene_id, start, end in intervals:
            fh.write(f"{gene_id}\t{start}\t{end}\n")


def _run_freebayes_dispatcher(
    cram: Path,
    ref_fa: Path,
    bed: Path,
    output_vcf: Path,
    n_threads: int = 4,
) -> None:
    """Run the Rust freebayes dispatcher with the current batch-size env."""
    if meteor_core is None:
        pytest.skip("meteor_core not installed")
    options = meteor_core.FreebayesOptions(
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
    )
    meteor_core.call_variants_parallel(
        str(cram),
        str(ref_fa),
        str(bed),
        FREEBAYES,
        options,
        n_threads=n_threads,
        output_path=str(output_vcf),
    )


def _normalize_vcf(input_vcf: Path, ref_fa: Path, output_vcf: Path) -> None:
    """Normalize and sort a VCF so records can be compared."""
    norm = subprocess.run(
        [BCFTOOLS, "norm", "-f", str(ref_fa), str(input_vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if norm.returncode != 0:
        raise RuntimeError(
            f"bcftools norm failed (exit {norm.returncode}): "
            f"{norm.stderr.decode('utf-8', errors='replace')}"
        )
    sort = subprocess.run(
        [
            BCFTOOLS,
            "sort",
            "-Oz",
            "-o",
            str(output_vcf),
            "-T",
            str(output_vcf.parent),
            "-",
        ],
        input=norm.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if sort.returncode != 0:
        raise RuntimeError(
            f"bcftools sort failed (exit {sort.returncode}): "
            f"{sort.stderr.decode('utf-8', errors='replace')}"
        )
    subprocess.run([TABIX, "-p", "vcf", str(output_vcf)], check=True)


def _vcf_records(path: Path) -> list[tuple[str, int, str, tuple[str, ...], Any]]:
    """Return a stable, comparable representation of VCF records."""
    with pysam.VariantFile(str(path)) as vcf:
        return [
            (
                rec.chrom,
                rec.pos,
                rec.ref,
                tuple(rec.alts) if rec.alts else (),
                rec.qual,
            )
            for rec in vcf
        ]


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #


def test_counter_tsv_byte_identity(tmp_path: Path) -> None:
    """Python and Rust counter paths must produce byte-identical TSV files."""
    paths = _env_paths()
    if meteor_core is None:
        pytest.skip("meteor_core not installed")

    cram = paths["METEOR_BENCH_CRAM"]
    mapping_dir = cram.parent
    ref_dir = paths["METEOR_BENCH_REF"]
    stage1_json, stage1_data, ref_json = _stage1_and_ref(mapping_dir, ref_dir)
    sample_name = stage1_json.name.replace("_census_stage_1.json", "")
    raw_cram = cram
    identity_threshold = stage1_data.get("counting", {}).get("identity_threshold", 0.95)

    py_count = tmp_path / "python.tsv.xz"
    rust_count = tmp_path / "rust.tsv.xz"
    py_stage1 = tmp_path / "py_stage1.json"
    rust_stage1 = tmp_path / "rust_stage1.json"
    py_stage1.write_text(json.dumps(stage1_data), encoding="utf-8")
    rust_stage1.write_text(json.dumps(stage1_data), encoding="utf-8")

    for use_rust, count_file, stage1_copy in (
        (False, py_count, py_stage1),
        (True, rust_count, rust_stage1),
    ):
        tmp = tmp_path / ("rust" if use_rust else "python")
        tmp.mkdir()
        meteor = Component(
            threads=8,
            fastq_dir=mapping_dir,
            mapping_dir=mapping_dir,
            ref_dir=ref_dir,
            tmp_path=tmp,
        )
        meteor.use_rust_counter = use_rust
        counter = Counter(
            meteor=meteor,
            counting_type="smart_shared",
            mapping_type="end-to-end",
            trim=80,
            identity_user=None,
            alignment_number=100,
            core_size=100,
            keep_all_alignments=False,
            keep_filtered_alignments=False,
            json_data={},
            identity_threshold=identity_threshold,
        )
        counter.launch_counting(
            raw_cram,
            mapping_dir / f"{sample_name}.cram",
            count_file,
            ref_json,
            json.loads(stage1_copy.read_text(encoding="utf-8")),
            stage1_copy,
        )

    with lzma.open(py_count, "rb") as fh:
        py_bytes = fh.read()
    with lzma.open(rust_count, "rb") as fh:
        rust_bytes = fh.read()

    if py_bytes != rust_bytes:
        diff_path = _evidence_dir() / f"{sample_name}_counter.tsv.diff"
        _write_text_diff(
            diff_path,
            "python.tsv.xz",
            py_bytes.decode("utf-8"),
            "rust.tsv.xz",
            rust_bytes.decode("utf-8"),
        )
        pytest.fail(f"counter TSV bytes differ; diff written to {diff_path}")


def test_vcf_record_equality_after_normalization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Rust freebayes dispatcher with batch_size=1 vs default must yield equal records."""
    paths = _env_paths()
    if meteor_core is None:
        pytest.skip("meteor_core not installed")
    if not shutil.which(FREEBAYES) or not shutil.which(BCFTOOLS):
        pytest.skip("freebayes and/or bcftools not available")

    cram = paths["METEOR_BENCH_CRAM"]
    ref_fa = _reference_fasta(paths["METEOR_BENCH_REF"], tmp_path)
    gene_limit = int(os.environ.get("METEOR_PARITY_VCF_GENES", "64"))
    intervals = _read_bed_intervals(paths["METEOR_BENCH_CATALOGUE"], limit=gene_limit)
    bed = tmp_path / "vcf_regions.bed"
    _write_bed(bed, intervals)

    default_vcf = tmp_path / "default.vcf"
    serial_vcf = tmp_path / "serial.vcf"

    monkeypatch.delenv("METEOR_FREEBAYES_BATCH_SIZE", raising=False)
    _run_freebayes_dispatcher(cram, ref_fa, bed, default_vcf, n_threads=2)

    monkeypatch.setenv("METEOR_FREEBAYES_BATCH_SIZE", "1")
    _run_freebayes_dispatcher(cram, ref_fa, bed, serial_vcf, n_threads=2)

    default_norm = tmp_path / "default.norm.vcf.gz"
    serial_norm = tmp_path / "serial.norm.vcf.gz"
    _normalize_vcf(default_vcf, ref_fa, default_norm)
    _normalize_vcf(serial_vcf, ref_fa, serial_norm)

    default_records = _vcf_records(default_norm)
    serial_records = _vcf_records(serial_norm)

    if default_records != serial_records:
        diff_path = _evidence_dir() / "vcf_serial_vs_default.diff"
        default_txt = "\n".join(repr(r) for r in default_records)
        serial_txt = "\n".join(repr(r) for r in serial_records)
        _write_text_diff(
            diff_path,
            "default_batch.norm.vcf.gz",
            default_txt,
            "batch_1.norm.vcf.gz",
            serial_txt,
        )
        pytest.fail(f"VCF records differ after norm|sort; diff written to {diff_path}")


def test_consensus_fasta_identity(tmp_path: Path) -> None:
    """Python and Rust create_consensus must produce byte-identical FASTA files."""
    paths = _env_paths()
    if meteor_core is None:
        pytest.skip("meteor_core not installed")
    if not shutil.which(FREEBAYES):
        pytest.skip("freebayes not available")

    cram = paths["METEOR_BENCH_CRAM"]
    ref_fa = _reference_fasta(paths["METEOR_BENCH_REF"], tmp_path)
    gene_limit = int(os.environ.get("METEOR_PARITY_CONSENSUS_GENES", "64"))
    intervals = _read_bed_intervals(paths["METEOR_BENCH_CATALOGUE"], limit=gene_limit)
    bed = tmp_path / "consensus_regions.bed"
    _write_bed(bed, intervals)

    vcf = tmp_path / "consensus.vcf"
    _run_freebayes_dispatcher(cram, ref_fa, bed, vcf, n_threads=2)

    sorted_vcf = tmp_path / "consensus.sorted.vcf.gz"
    subprocess.run(
        [
            BCFTOOLS,
            "sort",
            "-T",
            str(tmp_path),
            "-Oz",
            "-o",
            str(sorted_vcf),
            str(vcf),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    subprocess.run([TABIX, "-p", "vcf", str(sorted_vcf)], check=True)
    vcf = sorted_vcf

    matrix_rows = [
        {"gene_id": gene_id, "gene_length": end - start, "coverage": 1000}
        for gene_id, start, end in intervals
    ]
    matrix = pd.DataFrame(matrix_rows)
    matrix_file = tmp_path / "matrix.tsv.xz"
    matrix.to_csv(
        matrix_file,
        sep="\t",
        index=False,
        header=True,
        compression="xz",
    )

    meteor_rust = Component(
        threads=8, ref_dir=paths["METEOR_BENCH_REF"], tmp_path=tmp_path
    )
    meteor_rust.use_rust_variant_calling = True
    vc_rust = VariantCalling(
        meteor=meteor_rust,
        census={},
        max_depth=8000,
        min_depth=3,
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
        core_size=100,
    )
    vc_rust.matrix_file = matrix_file
    low_cov_sites, gene_ignore = vc_rust.filter_low_cov_sites(cram, ref_fa)

    py_consensus = tmp_path / "python_consensus.fasta.xz"
    rust_consensus = tmp_path / "rust_consensus.fasta.xz"

    meteor_py = Component(
        threads=8, ref_dir=paths["METEOR_BENCH_REF"], tmp_path=tmp_path
    )
    meteor_py.use_rust_variant_calling = False
    vc_py = VariantCalling(
        meteor=meteor_py,
        census={},
        max_depth=8000,
        min_depth=3,
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
        core_size=100,
    )
    vc_py.matrix_file = matrix_file
    vc_py.create_consensus(ref_fa, py_consensus, low_cov_sites, gene_ignore, vcf, bed)

    low_cov_list = (
        [
            (int(g), int(s), int(e))
            for g, s, e in low_cov_sites.reset_index()[
                ["gene_id", "startpos", "endpos"]
            ].itertuples(index=False, name=None)
        ]
        if not low_cov_sites.empty
        else []
    )
    ignore_list = (
        [
            (int(g), int(l))
            for g, l in gene_ignore.reset_index()[
                ["gene_id", "gene_length"]
            ].itertuples(index=False, name=None)
        ]
        if not gene_ignore.empty
        else []
    )
    records = meteor_core.create_consensus(
        str(ref_fa),
        str(vcf),
        str(bed),
        low_cov_list,
        ignore_list,
        0.1,
        Component.DEFAULT_GAP_CHAR,
    )
    with lzma.open(rust_consensus, "wt", preset=0) as fh:
        for gene_id, sequence in sorted(records, key=lambda item: item[0]):
            fh.write(f">{gene_id}\n{sequence}\n")

    with lzma.open(py_consensus, "rb") as fh:
        py_bytes = fh.read()
    with lzma.open(rust_consensus, "rb") as fh:
        rust_bytes = fh.read()

    if py_bytes != rust_bytes:
        diff_path = _evidence_dir() / "consensus_python_vs_rust.diff"
        _write_text_diff(
            diff_path,
            "python_consensus.fasta.xz",
            py_bytes.decode("utf-8"),
            "rust_consensus.fasta.xz",
            rust_bytes.decode("utf-8"),
        )
        pytest.fail(f"consensus FASTA bytes differ; diff written to {diff_path}")


def test_depth_array_equality(tmp_path: Path) -> None:
    """Python pysam pileup depths must equal Rust single-pass depths per gene."""
    paths = _env_paths()
    if meteor_core is None:
        pytest.skip("meteor_core not installed")

    cram = paths["METEOR_BENCH_CRAM"]
    ref_fa = _reference_fasta(paths["METEOR_BENCH_REF"], tmp_path)
    max_depth = 8000
    gene_limit = int(os.environ.get("METEOR_PARITY_DEPTH_GENES", "100"))
    intervals = _read_bed_intervals(paths["METEOR_BENCH_CATALOGUE"], limit=gene_limit)

    zero_based_intervals = [(gene_id, 0, end - start) for gene_id, start, end in intervals]
    rust_results = meteor_core.depth_per_gene(
        str(cram), str(ref_fa), zero_based_intervals, max_depth
    )

    vc = VariantCalling(
        meteor=Component(threads=1),
        census={},
        max_depth=max_depth,
        min_depth=3,
        min_snp_depth=1,
        min_frequency=0.1,
        ploidy=1,
        core_size=100,
    )

    failures: list[str] = []
    diff_dir = _evidence_dir()
    diff_dir.mkdir(parents=True, exist_ok=True)

    with pysam.AlignmentFile(
        str(cram), "rc", reference_filename=str(ref_fa)
    ) as cram_fh, pysam.FastaFile(str(ref_fa)) as fasta:
        for gene_depth in rust_results:
            gene_name = gene_depth.gene
            start, end = next((s, e) for g, s, e in intervals if g == gene_name)
            observed = list(gene_depth.depths)
            expected_dict = vc.count_reads_in_gene(
                cram_fh, gene_name, end - start, fasta
            )
            expected = [expected_dict.get(pos, 0) for pos in range(end - start)]
            if observed != expected:
                diff_path = diff_dir / f"{gene_name}.diff"
                with diff_path.open("w") as fh:
                    for pos, (obs, exp) in enumerate(zip(observed, expected)):
                        if obs != exp:
                            fh.write(f"{gene_name}\t{start + pos}\t{obs}\t{exp}\n")
                failures.append(gene_name)

    assert not failures, f"depth mismatch for genes: {failures[:10]}"
