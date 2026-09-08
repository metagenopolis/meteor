#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pandas",
#     "pysam",
#     "typer",
# ]
# ///

"""Benchmark Python vs Rust variant-calling hot paths using the eva71 test data."""

from __future__ import annotations

import json
import statistics
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter, process_time
from typing import NoReturn

import pandas as pd
import typer
from meteor.session import Component
from meteor.variantcalling import VariantCalling

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA_DIR = REPO_ROOT / "meteor" / "tests" / "test_variantcalling"
DEFAULT_REPORT = (
    REPO_ROOT
    / ".omo"
    / "evidence"
    / "meteor-rust-acceleration"
    / "benchmarks"
    / "hot_paths.json"
)


@dataclass
class RunMeasurement:
    wall_seconds: float
    cpu_seconds: float


@dataclass
class BenchEntry:
    wall_seconds_median: float
    cpu_seconds_median: float
    runs: list[RunMeasurement]


@dataclass
class Report:
    data_dir: str
    filter_low_cov_sites_python: BenchEntry
    filter_low_cov_sites_rust: BenchEntry
    create_consensus_python: BenchEntry
    create_consensus_rust: BenchEntry

    def to_dict(self) -> dict:
        return {
            "data_dir": self.data_dir,
            "filter_low_cov_sites_python": asdict(self.filter_low_cov_sites_python),
            "filter_low_cov_sites_rust": asdict(self.filter_low_cov_sites_rust),
            "create_consensus_python": asdict(self.create_consensus_python),
            "create_consensus_rust": asdict(self.create_consensus_rust),
        }


def _fatal(message: str) -> NoReturn:
    typer.echo(message, err=True)
    raise typer.Exit(1)


def _check_data(data_dir: Path) -> None:
    required = [
        "eva71/eva71_reference.json",
        "eva71/fasta/eva71.fasta.gz",
        "eva71/database/eva71.bed",
        "eva71_bench/eva71_bench.cram",
        "eva71_bench/eva71_bench.cram.crai",
        "eva71_bench/eva71_bench.tsv.xz",
        "eva71_bench/eva71_bench_census_stage_1.json",
        "eva71_bench/eva71_bench.vcf.gz",
    ]
    missing = [name for name in required if not (data_dir / name).is_file()]
    if missing:
        _fatal(f"Missing data file(s) in {data_dir}: {', '.join(missing)}")


def _build_vc(
    data_dir: Path, tmp_dir: Path, use_rust_variant_calling: bool
) -> VariantCalling:
    ref_json = json.loads(
        (data_dir / "eva71" / "eva71_reference.json").read_text(encoding="utf-8")
    )
    census_json = json.loads(
        (data_dir / "eva71_bench" / "eva71_bench_census_stage_1.json").read_text(
            encoding="utf-8"
        )
    )
    sample_info = census_json["sample_info"]
    stage3_dir = tmp_dir / sample_info["sample_name"]
    stage3_dir.mkdir(exist_ok=True, parents=True)

    meteor = Component
    meteor.ref_dir = data_dir / "eva71"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_dir
    meteor.tmp_path = tmp_dir
    meteor.strain_dir = tmp_dir / "out"
    meteor.use_rust_variant_calling = use_rust_variant_calling
    meteor.DEFAULT_GAP_CHAR = "?"

    data_dict = {
        "mapped_sample_dir": data_dir / "eva71_bench",
        "census": census_json,
        "directory": stage3_dir,
        "Stage3FileName": stage3_dir / "eva71_bench_census_stage_1.json",
        "reference": ref_json,
    }
    vc = VariantCalling(meteor, data_dict, 100, 3, 1, 0.01, 1, 100)
    vc.matrix_file = data_dir / "eva71_bench" / "eva71_bench.tsv.xz"
    return vc


def _median(runs: list[RunMeasurement]) -> BenchEntry:
    return BenchEntry(
        wall_seconds_median=statistics.median(r.wall_seconds for r in runs),
        cpu_seconds_median=statistics.median(r.cpu_seconds for r in runs),
        runs=runs,
    )


def _time_filter_low_cov_sites(
    data_dir: Path, tmp_dir: Path, use_rust_variant_calling: bool
) -> RunMeasurement:
    vc = _build_vc(data_dir, tmp_dir, use_rust_variant_calling)
    cram_file = data_dir / "eva71_bench" / "eva71_bench.cram"
    reference_file = data_dir / "eva71" / "fasta" / "eva71.fasta.gz"

    wall_start = perf_counter()
    cpu_start = process_time()
    vc.filter_low_cov_sites(cram_file, reference_file)
    wall_seconds = perf_counter() - wall_start
    cpu_seconds = process_time() - cpu_start
    return RunMeasurement(wall_seconds=wall_seconds, cpu_seconds=cpu_seconds)


def _time_create_consensus(
    data_dir: Path, tmp_dir: Path, use_rust_variant_calling: bool
) -> RunMeasurement:
    vc_py = _build_vc(data_dir, tmp_dir / "py", False)
    cram_file = data_dir / "eva71_bench" / "eva71_bench.cram"
    reference_file = data_dir / "eva71" / "fasta" / "eva71.fasta.gz"
    low_cov_sites, gene_ignore = vc_py.filter_low_cov_sites(cram_file, reference_file)

    vcf_file = data_dir / "eva71_bench" / "eva71_bench.vcf.gz"
    bed_file = data_dir / "eva71" / "database" / "eva71.bed"
    consensus_file = tmp_dir / "consensus.fasta.xz"

    vc = _build_vc(data_dir, tmp_dir, use_rust_variant_calling)
    wall_start = perf_counter()
    cpu_start = process_time()
    vc.create_consensus(
        reference_file,
        consensus_file,
        low_cov_sites,
        gene_ignore,
        vcf_file,
        bed_file,
    )
    wall_seconds = perf_counter() - wall_start
    cpu_seconds = process_time() - cpu_start
    return RunMeasurement(wall_seconds=wall_seconds, cpu_seconds=cpu_seconds)


@app.command()
def main(
    data_dir: Path = typer.Option(
        DEFAULT_DATA_DIR, "--data-dir", help="Directory containing eva71 test data."
    ),
    report: Path = typer.Option(
        DEFAULT_REPORT, "--report", help="Path to the JSON report file."
    ),
    runs: int = typer.Option(3, "--runs", help="Number of replicated timing runs."),
) -> None:
    """Time filter_low_cov_sites and create_consensus in Python and Rust."""
    if runs < 1:
        _fatal("--runs must be >= 1")

    data_dir = data_dir.resolve()
    _check_data(data_dir)

    tmp_dirs = [tempfile.mkdtemp(prefix="bench_vc_hot_") for _ in range(runs)]
    filter_python = _median(
        [_time_filter_low_cov_sites(data_dir, Path(tmp), False) for tmp in tmp_dirs]
    )
    tmp_dirs = [tempfile.mkdtemp(prefix="bench_vc_hot_") for _ in range(runs)]
    filter_rust = _median(
        [_time_filter_low_cov_sites(data_dir, Path(tmp), True) for tmp in tmp_dirs]
    )
    tmp_dirs = [tempfile.mkdtemp(prefix="bench_vc_hot_") for _ in range(runs)]
    consensus_python = _median(
        [_time_create_consensus(data_dir, Path(tmp), False) for tmp in tmp_dirs]
    )
    tmp_dirs = [tempfile.mkdtemp(prefix="bench_vc_hot_") for _ in range(runs)]
    consensus_rust = _median(
        [_time_create_consensus(data_dir, Path(tmp), True) for tmp in tmp_dirs]
    )

    report.parent.mkdir(parents=True, exist_ok=True)
    report_data = Report(
        data_dir=str(data_dir.resolve()),
        filter_low_cov_sites_python=filter_python,
        filter_low_cov_sites_rust=filter_rust,
        create_consensus_python=consensus_python,
        create_consensus_rust=consensus_rust,
    )
    report.write_text(json.dumps(report_data.to_dict(), indent=2), encoding="utf-8")

    typer.echo(f"Report written to {report}")
    typer.echo(
        f"filter_low_cov_sites: Python={filter_python.wall_seconds_median:.3f}s, "
        f"Rust={filter_rust.wall_seconds_median:.3f}s"
    )
    typer.echo(
        f"create_consensus: Python={consensus_python.wall_seconds_median:.3f}s, "
        f"Rust={consensus_rust.wall_seconds_median:.3f}s"
    )


if __name__ == "__main__":
    app()
