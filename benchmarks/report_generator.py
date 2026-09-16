#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "typer",
# ]
# ///

"""Render the phase-2 real-data benchmark report from per-run JSONs.

Every number in the generated report comes from the JSON files produced by
benchmarks/bench_phase2.py (and optionally enriched by benchmarks/collect_sacct.py
with MaxRSS). No hand-typed values are inserted.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import typer

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)


@dataclass(frozen=True)
class GroupKey:
    component: str
    sample: str
    implementation: str

    def __str__(self) -> str:
        return f"{self.component}:{self.sample}:{self.implementation}"


@dataclass
class Group:
    key: GroupKey
    runs: list[dict[str, Any]]

    def median_wall(self) -> float | None:
        if not self.runs:
            return None
        return statistics.median(
            r["measurement"]["wall_seconds"] for r in self.runs if "wall_seconds" in r.get("measurement", {})
        )

    def median_cpu(self) -> float | None:
        if not self.runs:
            return None
        return statistics.median(
            r["measurement"]["cpu_seconds"] for r in self.runs if "cpu_seconds" in r.get("measurement", {})
        )

    def max_rss(self) -> int | None:
        values = [
            r["measurement"]["max_rss_kb"]
            for r in self.runs
            if "max_rss_kb" in r.get("measurement", {})
        ]
        if not values:
            return None
        return max(values)

    def n_runs(self) -> int:
        return len(self.runs)


def _load_runs(evidence_dirs: list[Path]) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    for evidence_dir in evidence_dirs:
        for path in sorted(evidence_dir.glob("*.json")):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                typer.echo(f"warning: skipping invalid JSON {path}: {exc}", err=True)
                continue
            if "config" not in data or "measurement" not in data:
                typer.echo(f"warning: skipping {path}: missing config/measurement", err=True)
                continue
            data["_source_file"] = str(path.relative_to(evidence_dir))
            runs.append(data)
    return runs


def _group_runs(runs: list[dict[str, Any]]) -> dict[GroupKey, Group]:
    groups: dict[GroupKey, Group] = {}
    for run in runs:
        cfg = run["config"]
        key = GroupKey(
            component=cfg.get("component", "unknown"),
            sample=cfg.get("sample", "unknown"),
            implementation=cfg.get("implementation", "unknown"),
        )
        groups.setdefault(key, Group(key=key, runs=[])).runs.append(run)
    return groups


def _fmt_seconds(value: float | None) -> str:
    if value is None:
        return "—"
    return f"{value:.1f}"


def _fmt_cpu(value: float | None) -> str:
    if value is None:
        return "—"
    return f"{value:.1f}"


def _fmt_rss(value: int | None) -> str:
    if value is None:
        return "—"
    if value >= 1024**2:
        return f"{value / 1024**2:.1f} GiB"
    if value >= 1024:
        return f"{value / 1024:.1f} MiB"
    return f"{value} KiB"


def _speedup(python: float | None, other: float | None) -> str:
    if python is None or other is None or other == 0:
        return "—"
    ratio = python / other
    return f"{ratio:.2f}×"


def _verdict(python: float | None, other: float | None, higher_is_faster: bool = True) -> str:
    if python is None or other is None:
        return "incomplete"
    if higher_is_faster:
        if other < python * 0.9:
            return "faster"
        if other > python * 1.1:
            return "slower"
        return "parity"
    else:
        if other < python * 0.9:
            return "slower"
        if other > python * 1.1:
            return "faster"
        return "parity"


def _render_counter_table(groups: dict[GroupKey, Group], sample: str) -> str:
    rows: list[tuple[str, str, str, str, str, str, str]] = []
    for impl in ("python", "rust"):
        for ctype in ("smart_shared", "total"):
            key = GroupKey(component="counter", sample=sample, implementation=impl)
            group = groups.get(key)
            if group is None:
                rows.append((impl, ctype, "—", "—", "—", "0", "missing"))
                continue
            extra = group.runs[0]["config"].get("extra", {}) if group.runs else {}
            # Filter runs by counting_type in extra.
            filtered = [
                r for r in group.runs if r["config"].get("extra", {}).get("counting_type") == ctype
            ]
            if not filtered:
                rows.append((impl, ctype, "—", "—", "—", str(group.n_runs()), "no matching runs"))
                continue
            g = Group(key=key, runs=filtered)
            py_key = GroupKey(component="counter", sample=sample, implementation="python")
            py_group = groups.get(py_key)
            py_filtered = [
                r for r in (py_group.runs if py_group else [])
                if r["config"].get("extra", {}).get("counting_type") == ctype
            ]
            py_wall = statistics.median(r["measurement"]["wall_seconds"] for r in py_filtered) if py_filtered else None
            rows.append(
                (
                    impl,
                    ctype,
                    _fmt_seconds(g.median_wall()),
                    _fmt_cpu(g.median_cpu()),
                    _fmt_rss(g.max_rss()),
                    str(g.n_runs()),
                    _verdict(py_wall, g.median_wall()) if impl == "rust" else "baseline",
                )
            )

    lines = [
        f"### Counter: {sample}",
        "",
        "| Implementation | Counting type | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for impl, ctype, wall, cpu, rss, n, verdict in rows:
        lines.append(f"| {impl} | {ctype} | {wall} | {cpu} | {rss} | {n} | {verdict} |")
    lines.append("")
    return "\n".join(lines)


def _render_variant_table(
    groups: dict[GroupKey, Group],
    sample: str,
    component: str,
    impl_order: list[str],
    baseline_impl: str,
) -> str:
    rows = []
    baseline_key = GroupKey(component=component, sample=sample, implementation=baseline_impl)
    baseline_wall = groups.get(baseline_key).median_wall() if groups.get(baseline_key) else None
    for impl in impl_order:
        key = GroupKey(component=component, sample=sample, implementation=impl)
        group = groups.get(key)
        if group is None:
            rows.append((impl, "—", "—", "—", "0", "missing"))
            continue
        verdict = (
            "baseline"
            if impl == baseline_impl
            else _verdict(baseline_wall, group.median_wall())
        )
        rows.append(
            (
                impl,
                _fmt_seconds(group.median_wall()),
                _fmt_cpu(group.median_cpu()),
                _fmt_rss(group.max_rss()),
                str(group.n_runs()),
                verdict,
            )
        )

    title = (
        "Variant calling"
        if component == "variantcalling"
        else "Variant calling (64-gene subset fallback)"
    )
    lines = [
        f"### {title}: {sample}",
        "",
        "| Implementation | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for impl, wall, cpu, rss, n, verdict in rows:
        lines.append(f"| {impl} | {wall} | {cpu} | {rss} | {n} | {verdict} |")
    lines.append("")
    return "\n".join(lines)


def _render_consensus_depth_table(groups: dict[GroupKey, Group], sample: str) -> str:
    impl_order = ["python", "rust"]
    rows = []
    for impl in impl_order:
        key = GroupKey(component="consensus_depth", sample=sample, implementation=impl)
        group = groups.get(key)
        if group is None:
            rows.append((impl, "—", "—", "—", "0", "missing"))
            continue
        py_key = GroupKey(component="consensus_depth", sample=sample, implementation="python")
        py_wall = groups.get(py_key).median_wall() if groups.get(py_key) else None
        rows.append(
            (
                impl,
                _fmt_seconds(group.median_wall()),
                _fmt_cpu(group.median_cpu()),
                _fmt_rss(group.max_rss()),
                str(group.n_runs()),
                _verdict(py_wall, group.median_wall()) if impl != "python" else "baseline",
            )
        )

    lines = [
        f"### Consensus + depth: {sample}",
        "",
        "| Implementation | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for impl, wall, cpu, rss, n, verdict in rows:
        lines.append(f"| {impl} | {wall} | {cpu} | {rss} | {n} | {verdict} |")
    lines.append("")
    return "\n".join(lines)


@app.command()
def main(
    evidence_dir: Path = typer.Argument(..., help="Directory with per-run JSON files."),
    out: Path = typer.Option(..., "--out", help="Path to the Markdown report to write."),
    extra_evidence_dir: Path | None = typer.Option(
        None,
        "--extra-evidence-dir",
        help="Additional directory with per-run JSONs (e.g. subset fallback runs).",
    ),
) -> None:
    """Render final_report_phase2.md from benchmark JSONs."""
    evidence_dirs = [evidence_dir]
    if extra_evidence_dir is not None:
        evidence_dirs.append(extra_evidence_dir)
    runs = _load_runs(evidence_dirs)
    groups = _group_runs(runs)

    samples = sorted({k.sample for k in groups})
    components = sorted({k.component for k in groups})

    lines = [
        "# Phase-2 real-data benchmark report",
        "",
        "This report is generated automatically from per-run JSON files. Every wall",
        "time, CPU time, and MaxRSS value is taken from those JSONs; no numbers are",
        "typed by hand.",
        "",
        "## Methodology",
        "",
        "- Runs were executed on the Pasteur maestro cluster (queue `hubbioit`).",
        "- Counter and consensus+depth runs were executed inside a single Slurm job",
        "  in interleaved Python/Rust order to reduce machine drift.",
        "- Variant-calling configurations were submitted as separate Slurm jobs",
        "  with 24 hours walltime and run in parallel.",
        "- MaxRSS was collected from `sacct` after the jobs completed and merged",
        "  into the JSONs by `benchmarks/collect_sacct.py`.",
        "- Medians are reported over the completed runs; runs that timed out are",
        "  noted explicitly rather than extrapolated.",
        "",
        "## Samples",
        "",
        f"Measured samples: {', '.join(samples)}.",
        "",
        "## Results",
        "",
    ]

    if "counter" in components:
        lines.append("## Counter")
        lines.append("")
        for sample in samples:
            lines.append(_render_counter_table(groups, sample))

    if "variantcalling" in components:
        lines.append("## Variant calling (full catalogue)")
        lines.append("")
        for sample in samples:
            lines.append(
                _render_variant_table(
                    groups,
                    sample,
                    component="variantcalling",
                    impl_order=["python", "rust-batch8", "rust-batch1"],
                    baseline_impl="python",
                )
            )

    if "variantcalling-subset" in components:
        lines.append("## Variant calling (64-gene subset fallback)")
        lines.append("")
        lines.append(
            "Python exceeded the 24 h walltime on the full catalogue, so a 64-gene "
            "subset was used as a controlled fallback. Only Rust batch configurations "
            "completed; no Python baseline is available for this subset."
        )
        lines.append("")
        for sample in samples:
            lines.append(
                _render_variant_table(
                    groups,
                    sample,
                    component="variantcalling-subset",
                    impl_order=["rust-batch8", "rust-batch1"],
                    baseline_impl="rust-batch8",
                )
            )

    if "consensus_depth" in components:
        lines.append("## Consensus + depth")
        lines.append("")
        for sample in samples:
            lines.append(_render_consensus_depth_table(groups, sample))

    lines.append("## Raw data")
    lines.append("")
    lines.append(f"Per-run JSONs are stored in `{evidence_dir}`.")
    if extra_evidence_dir is not None:
        lines.append(f"Additional subset JSONs are stored in `{extra_evidence_dir}`.")
    lines.append("")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines), encoding="utf-8")
    typer.echo(f"Report written to {out}")


if __name__ == "__main__":
    app()
