# Rust acceleration benchmarks

These benchmarks compare the original Python implementation with the optional
Rust-accelerated paths on the repository test fixtures. Each configuration was
run 3 times and the median wall and CPU time is reported.

## Environment

- Machine: Apple M-series, macOS (local development machine)
- Python: 3.12
- Rust toolchain: stable
- Fixtures: `tests/data/fixtures` (120 genes, 10 MSPs, ~22 000 reads)

## Results

| Step | Implementation | Median wall (s) | Median CPU (s) | Notes |
|------------------|----------------|----------------:|---------------:|--------------------------------------------|
| `counter` | Python | 0.172 | 0.081 | Baseline Python hot loop |
| `counter` | Rust | 0.108 | 0.108 | ~1.6x wall-time speedup on this fixture |
| `variantcalling` | Python | 2.867 | 1.076 | Python freebayes dispatcher + consensus |
| `variantcalling` | Rust | 5.414 | 4.064 | Slower on the small fixture (see note) |

### Honest note on variant calling

The Rust-accelerated variant-calling path is **slower** than the Python path on
this tiny fixture. The Rust dispatcher adds serialisation overhead and the
fixture spends most of its time launching short `freebayes` processes, so the
fixed overhead dominates. The Rust helpers are expected to become beneficial
on larger catalogues with more regions and more reads per region, but this has
not been benchmarked yet.

## Raw per-run data

```json
{
  "meta": {
    "timestamp": "2026-09-08T08:05:16.475243+00:00",
    "fixture_dir": "/Users/aghozlan/workspace/meteor-rust-accel/tests/data/fixtures",
    "fixture_sizes": {
      "sample.cram": 1147931,
      "sample.cram.crai": 804,
      "reference.fa.gz": 75665,
      "reference.fa.gz.fai": 2471,
      "reference.json": 511,
      "annotation.tsv": 1374,
      "msp_map.tsv": 2465,
      "sample.tsv.xz": 784,
      "sample_census_stage_1.json": 659
    },
    "record_counts": {
      "n_genes": 120,
      "n_msps": 10,
      "n_reads": 22000
    }
  },
  "counter": {
    "wall_seconds_median": 0.10845762502867728,
    "cpu_seconds_median": 0.10836500000000004,
    "runs": [
      { "wall_seconds": 0.11395387497032061, "cpu_seconds": 0.10983100000000001 },
      { "wall_seconds": 0.10716558300191537, "cpu_seconds": 0.10714300000000004 },
      { "wall_seconds": 0.10845762502867728, "cpu_seconds": 0.10836500000000004 }
    ]
  },
  "variantcalling": {
    "wall_seconds_median": 5.414290707965847,
    "cpu_seconds_median": 4.064042000000001,
    "runs": [
      { "wall_seconds": 5.4513325419975445, "cpu_seconds": 4.0791450000000005 },
      { "wall_seconds": 5.377010792028159, "cpu_seconds": 4.033555000000001 },
      { "wall_seconds": 5.414290707965847, "cpu_seconds": 4.064042000000001 }
    ]
  }
}
```

## How to reproduce

With the Rust extension already built (`maturin develop -m rust/Cargo.toml`):

```bash
# Counter
python benchmarks/bench_counter.py --runs 3
python benchmarks/bench_counter.py --use-rust --runs 3

# Variant calling
python benchmarks/bench_variantcalling.py --runs 3
python benchmarks/bench_variantcalling.py --use-rust --runs 3
```

The scripts write their results to
`.omo/evidence/meteor-rust-acceleration/benchmarks/baseline.json`.
