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
import os
import subprocess
import sys
from pathlib import Path
from typing import NoReturn


def _sample_name_from_cram(cram: Path) -> str:
    """Derive the sample name from the CRAM filename."""
    return cram.stem

import pysam
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
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-acceleration" / "benchmarks"
)
REAL_DATA_EVIDENCE_DIR = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-speedup-phase2" / "benchmarks"
)

REAL_DATA_ENV_VARS: dict[str, str] = {
    "cram": "METEOR_BENCH_CRAM",
    "ref": "METEOR_BENCH_REF",
    "msp_map": "METEOR_BENCH_MSP_MAP",
    "catalogue": "METEOR_BENCH_CATALOGUE",
}


def _fatal(message: str) -> NoReturn:
    typer.echo(message, err=True)
    raise typer.Exit(1)


def _write_json(path: Path, meta, counter_runs, vc_runs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = Baseline(
        meta=meta,
        counter=_median_measurement(counter_runs),
        variantcalling=_median_measurement(vc_runs),
    )
    path.write_text(json.dumps(data.to_dict(), indent=2), encoding="utf-8")


def _real_data_mode_selected(
    cram: Path | None,
    ref: Path | None,
    msp_map: Path | None,
    catalogue: Path | None,
) -> bool:
    return any(arg is not None for arg in (cram, ref, msp_map, catalogue))


def _validate_real_data_env(
    cram: Path | None,
    ref: Path | None,
    msp_map: Path | None,
    catalogue: Path | None,
) -> dict[str, str]:
    """Return env-var values for the selected real-data args; fatal if any are missing."""
    selected = {
        name: arg
        for name, arg in (
            ("cram", cram),
            ("ref", ref),
            ("msp_map", msp_map),
            ("catalogue", catalogue),
        )
        if arg is not None
    }

    missing_env: list[str] = []
    values: dict[str, str] = {}
    for name in selected:
        env_name = REAL_DATA_ENV_VARS[name]
        value = os.environ.get(env_name)
        if not value:
            missing_env.append(env_name)
        else:
            values[name] = value

    if missing_env:
        _fatal(
            "Real-data benchmark mode requested but the following environment "
            f"variable(s) are missing or empty: {', '.join(missing_env)}. "
            "Set METEOR_BENCH_CRAM, METEOR_BENCH_REF, METEOR_BENCH_MSP_MAP "
            "(and METEOR_BENCH_CATALOGUE if --catalogue is used) and retry."
        )

    missing_files: list[str] = []
    for name, path_str in values.items():
        path = Path(path_str)
        if not path.is_file():
            missing_files.append(f"{REAL_DATA_ENV_VARS[name]}={path_str}")
        elif name == "cram":
            index_path = path.with_suffix(path.suffix + ".crai")
            if path.suffix == ".cram" and not index_path.is_file():
                missing_files.append(f"CRAM index {index_path}")
            try:
                with pysam.AlignmentFile(str(path), "rc") as _:
                    pass
            except Exception as exc:
                _fatal(
                    f"METEOR_BENCH_CRAM does not point to a readable CRAM file: {path} "
                    f"({exc})"
                )

    if missing_files:
        _fatal(
            "Real-data benchmark mode requested but the following file(s) are "
            f"missing or invalid: {', '.join(missing_files)}"
        )

    return values


def _run_fixture_mode(fixture_dir: Path, runs: int) -> None:
    """Run the original phase-1 interleaved fixture benchmark."""
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


def _run_real_data_mode(
    cram: Path,
    ref: Path,
    msp_map: Path,
    catalogue: Path | None,
    runs: int,
) -> None:
    """Validate real-data env vars and run the fast counter benchmark.

    The heavy lifting for real-data benchmarking lives in
    benchmarks/bench_phase2.py, which is designed to be driven from Slurm jobs
    for variant-calling and consensus+depth. This entry point runs the fast
    interleaved counter benchmark locally and points the user at the Slurm
    scripts for the slower components.
    """
    values = _validate_real_data_env(cram, ref, msp_map, catalogue)
    REAL_DATA_EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    marker = REAL_DATA_EVIDENCE_DIR / "realdata_mode_selected.marker"
    marker.write_text(
        f"Real-data mode selected with: {values}\nRuns requested: {runs}\n",
        encoding="utf-8",
    )

    sample_name = _sample_name_from_cram(Path(values["cram"]))
    mapping_dir = Path(values["cram"]).parent
    out_dir = REAL_DATA_EVIDENCE_DIR / "per_run"

    typer.echo("Running interleaved real-data counter benchmark ...")
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "benchmarks" / "bench_phase2.py"),
            "interleaved",
            str(mapping_dir),
            str(values["ref"]),
            sample_name,
            "--out-dir",
            str(out_dir),
            "--component",
            "counter",
            "--runs",
            str(runs),
        ],
        check=True,
    )

    typer.echo(
        "Counter benchmark complete. For variant-calling and consensus+depth "
        "benchmarks, submit the Slurm scripts generated in "
        f"{REAL_DATA_EVIDENCE_DIR / 'slurm'}."
    )


@app.command()
def main(
    fixture_dir: Path = typer.Option(
        DEFAULT_FIXTURE_DIR, "--fixture-dir", help="Directory containing fixture files."
    ),
    runs: int = typer.Option(5, "--runs", help="Number of interleaved runs per mode."),
    cram: Path | None = typer.Option(
        None,
        "--cram",
        help="Mode selector: benchmark real data (CRAM path from METEOR_BENCH_CRAM).",
    ),
    ref: Path | None = typer.Option(
        None,
        "--ref",
        help="Mode selector: benchmark real data (reference from METEOR_BENCH_REF).",
    ),
    msp_map: Path | None = typer.Option(
        None,
        "--msp-map",
        help="Mode selector: benchmark real data (MSP map from METEOR_BENCH_MSP_MAP).",
    ),
    catalogue: Path | None = typer.Option(
        None,
        "--catalogue",
        help="Mode selector: benchmark real data (catalogue from METEOR_BENCH_CATALOGUE).",
    ),
) -> None:
    """Run counter and variant-calling benchmarks in interleaved Python/Rust order.

    With no --cram/--ref/--msp-map/--catalogue arguments the fixture benchmark is
    run exactly as in phase 1. Passing any of those arguments selects real-data
    mode; the actual file paths are read from the METEOR_BENCH_* environment
    variables, which must be set.
    """
    if runs < 1:
        raise typer.Exit("--runs must be >= 1")

    if _real_data_mode_selected(cram, ref, msp_map, catalogue):
        _run_real_data_mode(cram, ref, msp_map, catalogue, runs)
    else:
        _run_fixture_mode(fixture_dir, runs)


if __name__ == "__main__":
    app()
