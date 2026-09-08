#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pysam",
#     "typer",
# ]
# ///

# ─── How to run ───
# 1. Install uv (if not installed):
#      curl -LsSf https://astral.sh/uv/install.sh | sh
# 2. Run directly:
#      uv run benchmarks/bench_variantcalling.py
# 3. With custom fixture or baseline paths:
#      uv run benchmarks/bench_variantcalling.py --fixture-dir ./fixtures --baseline ./baseline.json
# ──────────────────

"""Benchmark meteor's Python variant-calling path."""

from __future__ import annotations

import json
import statistics
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter, process_time
from typing import NoReturn

import typer
from meteor.session import Component
from meteor.variantcalling import VariantCalling

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "fixtures"
DEFAULT_BASELINE = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-acceleration" / "benchmarks" / "baseline.json"
)
REQUIRED_FIXTURES = [
    "sample.cram",
    "sample.cram.crai",
    "reference.fa.gz",
    "reference.fa.gz.fai",
    "reference.json",
    "annotation.tsv",
    "msp_map.tsv",
    "sample.tsv.xz",
    "sample_census_stage_1.json",
]


@dataclass
class RunMeasurement:
    """One timed execution of the variant-calling path."""

    wall_seconds: float
    cpu_seconds: float
    timestamp: str


@dataclass
class BenchEntry:
    """Aggregated benchmark entry for a single meteor module."""

    wall_seconds_median: float
    cpu_seconds_median: float
    runs: list[RunMeasurement]


@dataclass
class FixtureMeta:
    """Metadata describing the fixture data used for the benchmark."""

    timestamp: str
    fixture_dir: str
    fixture_sizes: dict[str, int]
    record_counts: dict[str, int]


@dataclass
class Baseline:
    """Root object written to the baseline evidence file."""

    meta: FixtureMeta
    counter: BenchEntry | None
    variantcalling: BenchEntry | None

    def to_dict(self) -> dict:
        """Convert to a JSON-serialisable dictionary."""
        return {
            "meta": asdict(self.meta),
            "counter": asdict(self.counter) if self.counter else None,
            "variantcalling": asdict(self.variantcalling) if self.variantcalling else None,
        }


def _fatal(message: str) -> NoReturn:
    typer.echo(message, err=True)
    raise typer.Exit(1)


def _check_fixtures(fixture_dir: Path) -> None:
    missing = [name for name in REQUIRED_FIXTURES if not (fixture_dir / name).is_file()]
    if missing:
        _fatal(
            f"Missing fixture file(s) in {fixture_dir}: {', '.join(missing)}. "
            "Run 'python tests/data/fixtures/make_fixtures.py' first."
        )


def _load_fixture_counts(fixture_dir: Path) -> dict[str, int]:
    info_path = fixture_dir / "fixture_info.json"
    if info_path.is_file():
        return dict(json.loads(info_path.read_text(encoding="utf-8")))
    return {}


def _bench_entry_from_dict(data: dict) -> BenchEntry:
    return BenchEntry(
        wall_seconds_median=data["wall_seconds_median"],
        cpu_seconds_median=data["cpu_seconds_median"],
        runs=[RunMeasurement(**run) for run in data["runs"]],
    )


def _load_baseline(path: Path) -> Baseline:
    if path.is_file():
        data = json.loads(path.read_text(encoding="utf-8"))
        return Baseline(
            meta=FixtureMeta(**data.get("meta", {})),
            counter=_bench_entry_from_dict(data["counter"]) if data.get("counter") else None,
            variantcalling=_bench_entry_from_dict(data["variantcalling"])
            if data.get("variantcalling")
            else None,
        )
    return Baseline(
        meta=FixtureMeta(
            timestamp="",
            fixture_dir="",
            fixture_sizes={},
            record_counts={},
        ),
        counter=None,
        variantcalling=None,
    )


def _build_meta(fixture_dir: Path) -> FixtureMeta:
    fixture_sizes = {
        name: (fixture_dir / name).stat().st_size for name in REQUIRED_FIXTURES
    }
    record_counts = _load_fixture_counts(fixture_dir)
    return FixtureMeta(
        timestamp=datetime.now(timezone.utc).isoformat(),
        fixture_dir=str(fixture_dir.resolve()),
        fixture_sizes=fixture_sizes,
        record_counts=record_counts,
    )


def _median_measurement(measurements: list[RunMeasurement]) -> BenchEntry:
    return BenchEntry(
        wall_seconds_median=statistics.median(m.wall_seconds for m in measurements),
        cpu_seconds_median=statistics.median(m.cpu_seconds for m in measurements),
        runs=measurements,
    )


def _run_variantcalling(fixture_dir: Path, run_index: int, use_rust: bool) -> RunMeasurement:
    ref_json = json.loads((fixture_dir / "reference.json").read_text(encoding="utf-8"))
    stage1_data = json.loads(
        (fixture_dir / "sample_census_stage_1.json").read_text(encoding="utf-8")
    )

    with tempfile.TemporaryDirectory(
        dir=fixture_dir, prefix=f"bench_vc_{run_index}_"
    ) as tmp_raw:
        tmp_dir = Path(tmp_raw)
        out_dir = tmp_dir / "out"
        out_dir.mkdir(parents=True)

        meteor = Component
        meteor.threads = 1
        meteor.tmp_path = tmp_dir
        meteor.tmp_dir = tmp_dir
        meteor.ref_dir = fixture_dir
        meteor.strain_dir = out_dir
        meteor.use_rust_variant_calling = use_rust

        census = {
            "mapped_sample_dir": fixture_dir,
            "directory": out_dir,
            "Stage3FileName": out_dir / "sample_census_stage_3.json",
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

        wall_start = perf_counter()
        cpu_start = process_time()
        variant_caller.execute()
        wall_seconds = perf_counter() - wall_start
        cpu_seconds = process_time() - cpu_start

    return RunMeasurement(
        wall_seconds=wall_seconds,
        cpu_seconds=cpu_seconds,
        timestamp=datetime.now(timezone.utc).isoformat(),
    )


@app.command()
def main(
    fixture_dir: Path = typer.Option(
        DEFAULT_FIXTURE_DIR, "--fixture-dir", help="Directory containing fixture files."
    ),
    baseline: Path = typer.Option(
        DEFAULT_BASELINE, "--baseline", help="Path to the baseline JSON file."
    ),
    runs: int = typer.Option(3, "--runs", help="Number of replicated timing runs."),
    use_rust: bool = typer.Option(
        False, "--use-rust", help="Benchmark the Rust-accelerated variant-calling helpers."
    ),
) -> None:
    """Time meteor.variantcalling.VariantCalling.execute on the fixture data."""
    if runs < 1:
        _fatal("--runs must be >= 1")

    fixture_dir = fixture_dir.resolve()
    _check_fixtures(fixture_dir)

    measurements = [_run_variantcalling(fixture_dir, i, use_rust) for i in range(runs)]
    bench_entry = _median_measurement(measurements)

    baseline.parent.mkdir(parents=True, exist_ok=True)
    baseline_data = _load_baseline(baseline)
    baseline_data.meta = _build_meta(fixture_dir)
    baseline_data.variantcalling = bench_entry

    baseline.write_text(json.dumps(baseline_data.to_dict(), indent=2), encoding="utf-8")
    mode = "Rust" if use_rust else "Python"
    typer.echo(
        f"Variant-calling ({mode}) median over {runs} runs: "
        f"wall={bench_entry.wall_seconds_median:.3f}s, "
        f"cpu={bench_entry.cpu_seconds_median:.3f}s -> {baseline}"
    )


if __name__ == "__main__":
    app()
