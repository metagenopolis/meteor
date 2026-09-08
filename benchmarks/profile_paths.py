#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pysam",
#     "typer",
# ]
# ///

"""Profile meteor's counter and variant-calling hot paths and emit a report."""

from __future__ import annotations

import cProfile
import importlib.util
import json
import pstats
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "fixtures"
EVIDENCE_DIR = (
    REPO_ROOT / ".omo" / "evidence" / "meteor-rust-acceleration" / "profiling"
)


def _import_bench_module(name: str) -> Any:
    path = REPO_ROOT / "benchmarks" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _profile_call(
    label: str, fixture_dir: Path, prof_path: Path, runner: Any
) -> pstats.Stats:
    prof_path.parent.mkdir(parents=True, exist_ok=True)
    profiler = cProfile.Profile()
    profiler.enable()
    runner(fixture_dir, 0)
    profiler.disable()
    profiler.dump_stats(str(prof_path))
    return pstats.Stats(str(prof_path))


def _top_functions(
    prof_path: Path, sort_key: str, limit: int = 5
) -> list[dict[str, Any]]:
    stats_full = pstats.Stats(str(prof_path))
    stats_full.calc_callees()
    total = sum(t[2] for t in stats_full.stats.values()) if stats_full.stats else 0.0

    stripped = pstats.Stats(str(prof_path)).strip_dirs()
    stripped.sort_stats(sort_key)

    results: list[dict[str, Any]] = []
    for (full_file, line, func), (cc, nc, tt, ct, callers) in stats_full.stats.items():
        if total:
            pct = tt / total * 100 if sort_key == "tottime" else ct / total * 100
        else:
            pct = 0.0
        stripped_name = Path(full_file).name
        results.append(
            {
                "full_file": full_file,
                "file": stripped_name,
                "line": line,
                "function": func,
                "ncalls": nc,
                "tottime": tt,
                "cumtime": ct,
                "percent": pct,
            }
        )
    if sort_key == "tottime":
        results.sort(key=lambda x: x["tottime"], reverse=True)
    else:
        results.sort(key=lambda x: x["cumtime"], reverse=True)
    return results[:limit]


def _module_in_meteor(func: dict[str, Any]) -> bool:
    return "meteor/" in str(func["full_file"]).replace(
        "\\", "/"
    ) or "site-packages/meteor" in str(func["full_file"]).replace("\\", "/")


def _decision(
    label: str, top_total: list[dict[str, Any]], top_cum: list[dict[str, Any]]
) -> str:
    meteor_share = sum(f["percent"] for f in top_total if _module_in_meteor(f))
    if meteor_share >= 20.0:
        return (
            f"GO for Rust acceleration of {label}: meteor-internal functions account for "
            f"~{meteor_share:.1f}% of exclusive CPU time."
        )
    if meteor_share >= 5.0:
        return (
            f"CONDITIONAL GO for {label}: meteor-internal functions account for "
            f"~{meteor_share:.1f}% of exclusive CPU time; verify with larger input before "
            "committing full implementation."
        )
    return (
        f"NO-GO for Rust acceleration of {label}: meteor-internal Python overhead is only "
        f"~{meteor_share:.1f}% of exclusive CPU time; the bottleneck is likely I/O or "
        "external subprocesses."
    )


def _format_row(func: dict[str, Any]) -> str:
    return (
        f"| `{func['function']}` | `{Path(func['file']).name}:{func['line']}` | "
        f"{func['ncalls']:,} | {func['tottime']:.4f}s | {func['cumtime']:.4f}s | "
        f"{func['percent']:.1f}% |"
    )


def main() -> None:
    fixture_dir = DEFAULT_FIXTURE_DIR
    counter_mod = _import_bench_module("bench_counter")
    vc_mod = _import_bench_module("bench_variantcalling")

    counter_stats = _profile_call(
        "counter", fixture_dir, EVIDENCE_DIR / "counter.prof", counter_mod._run_counter
    )
    vc_stats = _profile_call(
        "variantcalling",
        fixture_dir,
        EVIDENCE_DIR / "variantcalling.prof",
        vc_mod._run_variantcalling,
    )

    counter_prof = EVIDENCE_DIR / "counter.prof"
    vc_prof = EVIDENCE_DIR / "variantcalling.prof"
    counter_total = _top_functions(counter_prof, "tottime")
    counter_cum = _top_functions(counter_prof, "cumtime")
    vc_total = _top_functions(vc_prof, "tottime")
    vc_cum = _top_functions(vc_prof, "cumtime")

    report = f"""# Meteor hot-path profiling report

Fixture: `{fixture_dir}`

## counter.py (CRAM iteration / MSP counting)

### Top 5 by exclusive CPU time (`tottime`)

| Function | Location | Calls | Tottime | Cumtime | % of profile |
| --- | --- | --- | --- | --- | --- |
"""
    report += "\n".join(_format_row(f) for f in counter_total)
    report += "\n\n### Top 5 by cumulative CPU time (`cumtime`)\n\n| Function | Location | Calls | Tottime | Cumtime | % of profile |\n| --- | --- | --- | --- | --- | --- |\n"
    report += "\n".join(_format_row(f) for f in counter_cum)
    report += (
        f"\n\n**Decision:** {_decision('counter.py', counter_total, counter_cum)}\n"
    )

    report += f"""\n## variantcalling.py (coverage / consensus / freebayes dispatch)

### Top 5 by exclusive CPU time (`tottime`)

| Function | Location | Calls | Tottime | Cumtime | % of profile |
| --- | --- | --- | --- | --- | --- |
"""
    report += "\n".join(_format_row(f) for f in vc_total)
    report += "\n\n### Top 5 by cumulative CPU time (`cumtime`)\n\n| Function | Location | Calls | Tottime | Cumtime | % of profile |\n| --- | --- | --- | --- | --- | --- |\n"
    report += "\n".join(_format_row(f) for f in vc_cum)
    report += f"\n\n**Decision:** {_decision('variantcalling.py', vc_total, vc_cum)}\n"

    report += "\n## Raw profiles\n\n"
    report += f"- counter: `{EVIDENCE_DIR / 'counter.prof'}`\n"
    report += f"- variantcalling: `{EVIDENCE_DIR / 'variantcalling.prof'}`\n"

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    (EVIDENCE_DIR / "report.md").write_text(report, encoding="utf-8")

    summary = {
        "counter": {
            "top_tottime": counter_total,
            "top_cumtime": counter_cum,
            "decision": _decision("counter.py", counter_total, counter_cum),
        },
        "variantcalling": {
            "top_tottime": vc_total,
            "top_cumtime": vc_cum,
            "decision": _decision("variantcalling.py", vc_total, vc_cum),
        },
    }
    (EVIDENCE_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )

    print(f"Profile report written to {EVIDENCE_DIR / 'report.md'}")


if __name__ == "__main__":
    main()
