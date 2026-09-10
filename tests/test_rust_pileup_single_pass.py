"""Single-pass pileup parity tests for meteor_core.depth_per_gene."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

meteor_core = pytest.importorskip("meteor_core")

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"
MAX_DEPTH = 8000

GENES = ["1", "2", "3", "10", "42", "100"]


def _count_reads_in_gene_arrays(genes: list[tuple[str, int]]) -> dict[str, list[int]]:
    """Reference depth arrays via the existing per-gene Rust helper."""
    result: dict[str, list[int]] = {}
    for gene_name, gene_length in genes:
        counts = meteor_core.count_reads_in_gene(
            str(CRAM), str(REF), gene_name, gene_length, MAX_DEPTH
        )
        array = [counts.get(pos, 0) for pos in range(gene_length)]
        result[gene_name] = array
    return result


def _gene_length(gene_name: str) -> int:
    from pysam import FastaFile

    with FastaFile(str(REF)) as fasta:
        return fasta.get_reference_length(gene_name)


def test_depth_per_gene_parity_all_fixture_genes() -> None:
    """depth_per_gene must match per-gene count_reads_in_gene for every position."""
    genes = [(name, _gene_length(name)) for name in GENES]
    expected = _count_reads_in_gene_arrays(genes)
    gene_intervals = [(name, 0, length) for name, length in genes]
    results = meteor_core.depth_per_gene(str(CRAM), str(REF), gene_intervals, MAX_DEPTH)

    assert len(results) == len(genes)
    for gene_depth in results:
        assert list(gene_depth.depths) == expected[gene_depth.gene]


def test_depth_per_gene_overlapping_intervals() -> None:
    """Overlapping intervals on the same gene accumulate independently."""
    gene_name = "1"
    gene_length = _gene_length(gene_name)
    half = gene_length // 2
    intervals = [
        (gene_name, 0, gene_length),
        (gene_name, 0, half),
        (gene_name, half, gene_length),
    ]
    results = meteor_core.depth_per_gene(str(CRAM), str(REF), intervals, MAX_DEPTH)

    assert len(results) == 3
    assert list(results[0].depths) == list(results[1].depths) + list(results[2].depths)


def test_depth_per_gene_uncovered_interval() -> None:
    """An interval with zero read coverage returns a zero-filled array of the right length."""
    gene_name = "1"
    gene_length = _gene_length(gene_name)
    counts = meteor_core.count_reads_in_gene(
        str(CRAM), str(REF), gene_name, gene_length, MAX_DEPTH
    )
    window_start = next(
        pos
        for pos in range(gene_length - 10)
        if all(counts.get(pos + offset, 0) == 0 for offset in range(10))
    )
    window_end = window_start + 10
    results = meteor_core.depth_per_gene(
        str(CRAM), str(REF), [(gene_name, window_start, window_end)], MAX_DEPTH
    )

    assert len(results) == 1
    assert len(results[0].depths) == window_end - window_start
    assert all(d == 0 for d in results[0].depths)


def test_depth_per_gene_max_depth_cap() -> None:
    """A high-coverage position is capped at max_depth in the single-pass path."""
    gene_name = "1"
    gene_length = _gene_length(gene_name)
    cap = 5
    full_depth = meteor_core.count_reads_in_gene(
        str(CRAM), str(REF), gene_name, gene_length, MAX_DEPTH
    )
    max_full = max(full_depth.values()) if full_depth else 0
    assert max_full > cap, "fixture max depth is below cap; test is vacuous"

    single_pass = meteor_core.depth_per_gene(
        str(CRAM), str(REF), [(gene_name, 0, gene_length)], cap
    )
    observed = list(single_pass[0].depths)
    assert all(d <= cap for d in observed)
    assert any(d == cap for d in observed)


def test_depth_per_gene_missing_cram() -> None:
    """A missing CRAM path raises a clear error."""
    missing = str(FIXTURE_DIR / "does_not_exist.cram")
    with pytest.raises(Exception):  # noqa: B017
        meteor_core.depth_per_gene(missing, str(REF), [("1", 0, 10)], MAX_DEPTH)


@pytest.mark.skipif(
    os.environ.get("METEOR_REAL_BENCH_DATA") != "1",
    reason="METEOR_REAL_BENCH_DATA != 1",
)
def test_depth_per_gene_real_data_parity() -> None:
    """Compare single-pass depths to per-gene depths on real benchmark data.

    Failures write per-gene diff files to
    .omo/evidence/meteor-rust-speedup-phase2/realdata/.
    """
    import csv

    cram = Path(os.environ["METEOR_BENCH_CRAM"])
    ref = Path(os.environ["METEOR_BENCH_REF"])
    catalogue = Path(os.environ["METEOR_BENCH_CATALOGUE"])
    if not cram.exists() or not ref.exists() or not catalogue.exists():
        pytest.skip("METEOR_BENCH_CRAM, METEOR_BENCH_REF, or METEOR_BENCH_CATALOGUE not found")

    # Read gene intervals from the catalogue BED (gene_id, start, end).
    intervals: list[tuple[str, int, int]] = []
    with catalogue.open() as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 3:
                continue
            intervals.append((row[0], int(row[1]), int(row[2])))

    results = meteor_core.depth_per_gene(str(cram), str(ref), intervals, MAX_DEPTH)

    diff_dir = (
        Path(__file__).resolve().parent.parent
        / ".omo"
        / "evidence"
        / "meteor-rust-speedup-phase2"
        / "realdata"
    )
    diff_dir.mkdir(parents=True, exist_ok=True)

    failures: list[str] = []
    for gene_depth in results:
        gene_name = gene_depth.gene
        start, end = next((s, e) for g, s, e in intervals if g == gene_name)
        expected = meteor_core.count_reads_in_gene(
            str(cram), str(ref), gene_name, end - start, MAX_DEPTH
        )
        expected_array = [expected.get(pos, 0) for pos in range(start, end)]
        observed_array = list(gene_depth.depths)
        if observed_array != expected_array:
            diff_path = diff_dir / f"{gene_name}.diff"
            with diff_path.open("w") as f:
                for pos, (obs, exp) in enumerate(zip(observed_array, expected_array)):
                    if obs != exp:
                        f.write(f"{gene_name}\t{start + pos}\t{obs}\t{exp}\n")
            failures.append(gene_name)

    assert not failures, f"depth mismatch for genes: {failures[:10]}..."
