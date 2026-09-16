#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer",
# ]
# ///

"""Merge MaxRSS from sacct into per-run benchmark JSONs.

Reads every benchmark JSON in the input directory, queries sacct for the
recorded SLURM_JOB_ID, and writes the MaxRSS value (in KiB) back into the JSON.
Runs that have no job_id or whose job is not found in sacct are left unchanged
but annotated with a note.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from typing import Any

import typer

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)


def _parse_max_rss(value: str | None) -> int | None:
    if not value or value == "-":
        return None
    match = re.match(r"([0-9.]+)([KMGTP]?)", value.strip())
    if not match:
        return None
    number = float(match.group(1))
    unit = match.group(2) or "K"
    multipliers = {"K": 1, "M": 1024, "G": 1024**2, "T": 1024**3, "P": 1024**4}
    return int(number * multipliers.get(unit, 1))


def _fetch_sacct(job_ids: set[str]) -> dict[str, int]:
    if not job_ids:
        return {}
    ids = ",".join(sorted(job_ids))
    result = subprocess.run(
        [
            "sacct",
            "--jobs",
            ids,
            "--format",
            "JobID,MaxRSS",
            "--noheader",
            "--parsable2",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    rss_by_job: dict[str, int] = {}
    if result.returncode != 0:
        return rss_by_job
    for line in result.stdout.splitlines():
        parts = line.split("|")
        if len(parts) < 2:
            continue
        full_job_id = parts[0].strip()
        job_id = full_job_id.split(".")[0]
        rss = _parse_max_rss(parts[1].strip())
        if rss is not None and job_id:
            rss_by_job[job_id] = max(rss_by_job.get(job_id, 0), rss)
    return rss_by_job


def _collect(in_dir: Path, dry_run: bool) -> None:
    json_paths = sorted(in_dir.glob("*.json"))
    if not json_paths:
        typer.echo(f"No JSON files found in {in_dir}")
        raise typer.Exit(1)

    job_ids: set[str] = set()
    for path in json_paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        job_id = data.get("measurement", {}).get("job_id")
        if job_id:
            job_ids.add(str(job_id))

    rss_by_job = _fetch_sacct(job_ids)

    for path in json_paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        measurement: dict[str, Any] = data.setdefault("measurement", {})
        job_id = measurement.get("job_id")
        notes: list[str] = data.setdefault("notes", [])
        if not job_id:
            if "no job_id" not in notes:
                notes.append("no job_id")
        else:
            rss = rss_by_job.get(str(job_id))
            if rss is not None:
                measurement["max_rss_kb"] = rss
            else:
                note = f"sacct returned no MaxRSS for {job_id}"
                if note not in notes:
                    notes.append(note)
        if not dry_run:
            path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        typer.echo(f"updated {path.name}")


@app.command()
def main(
    in_dir: Path = typer.Argument(..., help="Directory containing per-run JSON files."),
    dry_run: bool = typer.Option(False, "--dry-run"),
) -> None:
    """Merge sacct MaxRSS into benchmark JSONs."""
    _collect(in_dir, dry_run)


if __name__ == "__main__":
    app()
