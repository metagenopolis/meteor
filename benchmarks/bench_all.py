#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pysam",
#     "typer",
# ]
# ///

"""Interleaved benchmark driver for the Rust-accelerated meteor paths.

Runs Python and Rust variants in an alternating order so machine drift affects
both implementations equally. Writes `python.json` and `rust.json` in the
evidence directory, each containing 5 runs for `counter` and `variantcalling`.
"""

from __future__ import annotations

import json
from pathlib import Path

import typer

from bench_counter import (
    Baseline,
    RunMeasurement,
    _build_meta,
    _median_measurement,
    _run_counter,
)
from bench_variantcalling import _run_variantcalling

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "fixtures"
EVIDENCE_DIR = (
    REPO_ROOT
    / ".omo"
    / "evidence"
    / "meteor-rust-acceleration"
    / "benchmarks"
)


def _write_json(path: Path, meta, counter_runs, vc_runs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = Baseline(
        meta=meta,
        counter=_median_measurement(counter_runs),
        variantcalling=_median_measurement(vc_runs),
    )
    path.write_text(json.dumps(data.to_dict(), indent=2), encoding="utf-8")


@app.command()
def main(
    fixture_dir: Path = typer.Option(
        DEFAULT_FIXTURE_DIR, "--fixture-dir", help="Directory containing fixture files."
    ),
    runs: int = typer.Option(5, "--runs", help="Number of interleaved runs per mode."),
) -> None:
    """Run counter and variant-calling benchmarks in interleaved Python/Rust order."""
    if runs < 1:
        raise typer.Exit("--runs must be >= 1")

    fixture_dir = fixture_dir.resolve()
    meta = _build_meta(fixture_dir)

    python_counter_runs: list[RunMeasurement] = []
    rust_counter_runs: list[RunMeasurement] = []
    python_vc_runs: list[RunMeasurement] = []
    rust_vc_runs: list[RunMeasurement] = []

    for i in range(runs):
        typer.echo(f"Run {i + 1}/{runs} ...")
        python_counter_runs.append(_run_counter(fixture_dir, i, use_rust=False))
        rust_counter_runs.append(_run_counter(fixture_dir, i, use_rust=True))
        python_vc_runs.append(_run_variantcalling(fixture_dir, i, use_rust=False))
        rust_vc_runs.append(_run_variantcalling(fixture_dir, i, use_rust=True))

    _write_json(EVIDENCE_DIR / "python.json", meta, python_counter_runs, python_vc_runs)
    _write_json(EVIDENCE_DIR / "rust.json", meta, rust_counter_runs, rust_vc_runs)

    py_counter = _median_measurement(python_counter_runs)
    rust_counter = _median_measurement(rust_counter_runs)
    py_vc = _median_measurement(python_vc_runs)
    rust_vc = _median_measurement(rust_vc_runs)

    typer.echo(
        f"Counter  — Python: wall={py_counter.wall_seconds_median:.3f}s, "
        f"cpu={py_counter.cpu_seconds_median:.3f}s; "
        f"Rust: wall={rust_counter.wall_seconds_median:.3f}s, "
        f"cpu={rust_counter.cpu_seconds_median:.3f}s"
    )
    typer.echo(
        f"Variant  — Python: wall={py_vc.wall_seconds_median:.3f}s, "
        f"cpu={py_vc.cpu_seconds_median:.3f}s; "
        f"Rust: wall={rust_vc.wall_seconds_median:.3f}s, "
        f"cpu={rust_vc.cpu_seconds_median:.3f}s"
    )
    typer.echo(f"Results written to {EVIDENCE_DIR}/python.json and rust.json")


if __name__ == "__main__":
    app()
