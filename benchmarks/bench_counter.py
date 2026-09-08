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
#      uv run benchmarks/bench_counter.py
# 3. With custom fixture or baseline paths:
#      uv run benchmarks/bench_counter.py --fixture-dir ./fixtures --baseline ./baseline.json
# ──────────────────

"""Benchmark meteor's Python CRAM counting hot loop."""

from __future__ import annotations

import json
import shutil
import statistics
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter, process_time
from typing import NoReturn

import typer
from meteor.counter import Counter
from meteor.session import Component

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "fixtures"
DEFAULT_BASELINE = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-acceleration" / "benchmarks" / "baseline.json"
)
REQUIRED_FIXTURES = [
    "sample.cram",
    "sample.cram.crai",
    "reference.fa",
    "reference.fa.fai",
    "reference.json",
    "sample_census_stage_1.json",
]


@dataclass
class RunMeasurement:
    """One timed execution of the counter hot loop."""

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


def _run_counter(fixture_dir: Path, run_index: int, use_rust: bool) -> RunMeasurement:
    ref_json = json.loads((fixture_dir / "reference.json").read_text(encoding="utf-8"))
    stage1_data = json.loads(
        (fixture_dir / "sample_census_stage_1.json").read_text(encoding="utf-8")
    )

    with tempfile.TemporaryDirectory(
        dir=fixture_dir, prefix=f"bench_counter_{run_index}_"
    ) as tmp_raw:
        tmp_dir = Path(tmp_raw)
        meteor = Component
        meteor.threads = 1
        meteor.tmp_path = tmp_dir
        meteor.tmp_dir = tmp_dir
        meteor.mapping_dir = fixture_dir
        meteor.fastq_dir = fixture_dir
        meteor.ref_dir = fixture_dir
        meteor.use_rust_counter = use_rust

        counter = Counter(
            meteor,
            "smart_shared",
            "end-to-end",
            80,
            None,
            100,
            100,
        )

        wall_start = perf_counter()
        cpu_start = process_time()
        counter.launch_counting(
            fixture_dir / "sample.cram",
            tmp_dir / f"strain_{run_index}.cram",
            tmp_dir / f"sample_{run_index}.tsv.xz",
            ref_json,
            stage1_data,
            tmp_dir / f"sample_{run_index}_census_stage_1.json",
        )
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
    output: Path = typer.Option(
        None,
        "--output",
        help="Path to the JSON file to write. Defaults to python.json or rust.json under the evidence directory.",
    ),
    runs: int = typer.Option(3, "--runs", help="Number of replicated timing runs."),
    mode: str = typer.Option(
        "python",
        "--mode",
        help="Which implementation to benchmark.",
        case_sensitive=False,
    ),
) -> None:
    """Time meteor.counter.launch_counting on the fixture CRAM."""
    if runs < 1:
        _fatal("--runs must be >= 1")
    mode = mode.lower()
    if mode not in ("python", "rust"):
        _fatal("--mode must be 'python' or 'rust'")
    use_rust = mode == "rust"

    fixture_dir = fixture_dir.resolve()
    _check_fixtures(fixture_dir)

    if output is None:
        output = (
            Path(__file__).resolve().parent.parent
            / ".omo"
            / "evidence"
            / "meteor-rust-acceleration"
            / "benchmarks"
            / f"{mode}.json"
        )

    measurements = [_run_counter(fixture_dir, i, use_rust) for i in range(runs)]
    bench_entry = _median_measurement(measurements)

    output.parent.mkdir(parents=True, exist_ok=True)
    data = _load_baseline(output)
    data.meta = _build_meta(fixture_dir)
    data.counter = bench_entry

    output.write_text(json.dumps(data.to_dict(), indent=2), encoding="utf-8")
    mode_label = mode.capitalize()
    typer.echo(
        f"Counter ({mode_label}) median over {runs} runs: "
        f"wall={bench_entry.wall_seconds_median:.3f}s, "
        f"cpu={bench_entry.cpu_seconds_median:.3f}s -> {output}"
    )


if __name__ == "__main__":
    app()
