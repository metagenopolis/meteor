# Rust acceleration benchmarks

These benchmarks compare the original Python implementation with the optional
Rust-accelerated paths on the repository test fixtures. Each configuration was
run 3 times and the median wall and CPU time is reported.

The numbers below are derived directly from
`.omo/evidence/meteor-rust-acceleration/benchmarks/python.json` and
`rust.json` (see raw JSON blocks at the end of this file).

## Environment

- Machine: Apple M-series, macOS (local development machine)
- Python: 3.12
- Rust toolchain: stable
- Fixtures: `tests/data/fixtures` (120 genes, 10 MSPs, ~22 000 reads)

## Results

| Step | Implementation | Median wall (s) | Median CPU (s) | Notes |
|------|----------------|----------------:|---------------:|-------|
| `counter` | Python | 0.152 | 0.110 | Baseline Python hot loop |
| `counter` | Rust | 0.166 | 0.166 | Equivalent on this small fixture |
| `variantcalling` | Python | 8.236 | 2.976 | Python freebayes dispatcher + consensus |
| `variantcalling` | Rust | 15.465 | 10.481 | Slower on the small fixture (see note) |

### Honest notes

- `counter`: the Rust path is in the same ballpark as Python on this tiny fixture;
  the overhead of crossing the Python/Rust boundary masks any raw-loop speed-up.
- `variantcalling`: the Rust-accelerated path is **slower** than the Python path on
  this tiny fixture. The Rust dispatcher adds serialisation overhead and the
  fixture spends most of its time launching short `freebayes` processes, so the
  fixed overhead dominates. The Rust helpers are expected to become beneficial
  on larger catalogues with more regions and more reads per region, but this has
  not been benchmarked yet.

## Raw per-run data — Python

```json
{
  "counter": {
    "wall_seconds_median": 0.1520134579623118,
    "cpu_seconds_median": 0.109591,
    "runs": [
      { "wall_seconds": 0.17252983298385516, "cpu_seconds": 0.10992399999999997 },
      { "wall_seconds": 0.1510044580209069, "cpu_seconds": 0.109591 },
      { "wall_seconds": 0.1520134579623118, "cpu_seconds": 0.10910300000000006 }
    ]
  },
  "variantcalling": {
    "wall_seconds_median": 8.236228374997154,
    "cpu_seconds_median": 2.975575,
    "runs": [
      { "wall_seconds": 8.236228374997154, "cpu_seconds": 2.862864 },
      { "wall_seconds": 8.957224791985936, "cpu_seconds": 3.0513979999999994 },
      { "wall_seconds": 7.602691125008278, "cpu_seconds": 2.975575 }
    ]
  }
}
```

## Raw per-run data — Rust

```json
{
  "counter": {
    "wall_seconds_median": 0.16595620795851573,
    "cpu_seconds_median": 0.165535,
    "runs": [
      { "wall_seconds": 0.158267209015321, "cpu_seconds": 0.15305899999999995 },
      { "wall_seconds": 0.16595620795851573, "cpu_seconds": 0.165535 },
      { "wall_seconds": 0.1975272920099087, "cpu_seconds": 0.19591900000000007 }
    ]
  },
  "variantcalling": {
    "wall_seconds_median": 15.464665875013452,
    "cpu_seconds_median": 10.480587999999997,
    "runs": [
      { "wall_seconds": 19.79193554102676, "cpu_seconds": 12.292219 },
      { "wall_seconds": 14.777322875044774, "cpu_seconds": 10.378901 },
      { "wall_seconds": 15.464665875013452, "cpu_seconds": 10.480587999999997 }
    ]
  }
}
```

## How to reproduce

With the Rust extension already built (`maturin develop -m rust/Cargo.toml`):

```bash
# Counter
python benchmarks/bench_counter.py --mode python --runs 3
python benchmarks/bench_counter.py --mode rust --runs 3

# Variant calling
python benchmarks/bench_variantcalling.py --mode python --runs 3
python benchmarks/bench_variantcalling.py --mode rust --runs 3
```

The scripts write their results to
`.omo/evidence/meteor-rust-acceleration/benchmarks/python.json` and
`rust.json`.
