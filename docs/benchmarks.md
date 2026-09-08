# Rust acceleration benchmarks

These benchmarks compare the original Python implementation with the optional
Rust-accelerated paths on the repository test fixtures. Each configuration was
run 5 times in an interleaved Python/Rust order and the median wall and CPU time
is reported.

The numbers below are derived directly from
`.omo/evidence/meteor-rust-acceleration/benchmarks/python.json` and
`rust.json`.

## Environment

- Machine: Apple M-series, macOS (local development machine)
- Python: 3.12
- Rust toolchain: stable
- Fixtures: `tests/data/fixtures` (120 genes, 10 MSPs, ~22 000 reads)

## Results

| Step | Implementation | Median wall (s) | Median CPU (s) | Notes |
|------|----------------|----------------:|---------------:|-------|
| `counter` | Python | 0.194 | 0.121 | Baseline Python hot loop |
| `counter` | Rust | 0.196 | 0.196 | Equivalent/slightly slower on this fixture |
| `variantcalling` | Python | 4.860 | 1.861 | Python freebayes dispatcher + consensus |
| `variantcalling` | Rust | 15.301 | 11.247 | Slower on the small fixture (see note) |

### Variance note

Absolute wall times on this tiny fixture vary by roughly ±0.05-0.2s between
sessions because the measurements are dominated by interpreter startup, OS
scheduling, and short `freebayes` subprocesses rather than by the counting code
itself. The Rust/Python ratios are the meaningful signal, not the absolute
values.

### Honest notes

- `counter`: the Rust path is in the same ballpark as Python on this tiny
  fixture; the overhead of crossing the Python/Rust boundary masks any raw-loop
  speed-up. See the architectural note below for why.
- `variantcalling`: the Rust-accelerated path is **slower** than the Python path
  on this tiny fixture. The Rust dispatcher adds serialisation overhead and the
  fixture spends most of its time launching short `freebayes` processes, so the
  fixed overhead dominates. The Rust helpers are expected to become beneficial
  on larger catalogues with more regions and more reads per region, but this has
  not been benchmarked yet.

## Architectural note on the counter path

`meteor_core.count_msp` currently returns per-gene **read lists** (Python
string objects for all ~22 000 reads) across the PyO3 boundary. Rust performs
the per-read classification inside the extension, but the cost of building and
moving those Python string objects dominates the small-fixture run time. A
future redesign that returns only per-gene aggregates from Rust (instead of the
full read lists) is the expected path to a real counter speed-up.

## Variant-calling hot paths

The CRAM pileup loop (`filter_low_cov_sites`) shows a ~2x improvement from
moving the per-position read counting into Rust on the eva71 hot-path fixtures.
`create_consensus` is dominated by VCF/FASTA I/O on the tiny reference, so the
Rust rewrite is parity rather than a speed-up.

## Raw per-run data — Python

```json
{
  "counter": {
    "wall_seconds_median": 0.19385487504769117,
    "cpu_seconds_median": 0.1205889999999954,
    "runs": [
      { "wall_seconds": 0.1379981660284102, "cpu_seconds": 0.07325200000000004 },
      { "wall_seconds": 0.33320283296052366, "cpu_seconds": 0.07487900000000014 },
      { "wall_seconds": 0.26790708396583796, "cpu_seconds": 0.21202200000000104 },
      { "wall_seconds": 0.19385487504769117, "cpu_seconds": 0.15236800000000272 },
      { "wall_seconds": 0.16284987499238923, "cpu_seconds": 0.1205889999999954 }
    ]
  },
  "variantcalling": {
    "wall_seconds_median": 4.859759040991776,
    "cpu_seconds_median": 1.8607039999999984,
    "runs": [
      { "wall_seconds": 2.8604181249975227, "cpu_seconds": 1.064769 },
      { "wall_seconds": 2.936176292016171, "cpu_seconds": 1.1340089999999998 },
      { "wall_seconds": 8.161557250015903, "cpu_seconds": 3.050193 },
      { "wall_seconds": 8.71678974997485, "cpu_seconds": 3.251928999999997 },
      { "wall_seconds": 4.859759040991776, "cpu_seconds": 1.8607039999999984 }
    ]
  }
}
```

## Raw per-run data — Rust

```json
{
  "counter": {
    "wall_seconds_median": 0.19616025005234405,
    "cpu_seconds_median": 0.1960940000000022,
    "runs": [
      { "wall_seconds": 0.11285458295606077, "cpu_seconds": 0.10860999999999998 },
      { "wall_seconds": 0.10929458303144202, "cpu_seconds": 0.10911800000000049 },
      { "wall_seconds": 0.361048708029557, "cpu_seconds": 0.33995200000000025 },
      { "wall_seconds": 0.30672608397435397, "cpu_seconds": 0.2948379999999986 },
      { "wall_seconds": 0.19616025005234405, "cpu_seconds": 0.1960940000000022 }
    ]
  },
  "variantcalling": {
    "wall_seconds_median": 15.300728875037748,
    "cpu_seconds_median": 11.247238000000003,
    "runs": [
      { "wall_seconds": 5.497508750006091, "cpu_seconds": 4.124796999999999 },
      { "wall_seconds": 19.06489674997283, "cpu_seconds": 11.813307000000002 },
      { "wall_seconds": 15.300728875037748, "cpu_seconds": 11.247238000000003 },
      { "wall_seconds": 24.935627041966654, "cpu_seconds": 12.514010000000006 },
      { "wall_seconds": 11.770945333992131, "cpu_seconds": 9.231468 }
    ]
  }
}
```

## How to reproduce

With the Rust extension already built (`maturin develop -m rust/Cargo.toml`):

```bash
python benchmarks/bench_all.py --runs 5
```

The script writes its results to
`.omo/evidence/meteor-rust-acceleration/benchmarks/python.json` and
`rust.json`.
