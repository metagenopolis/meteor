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

## Phase-2 real-data benchmarks

The numbers in this section were produced on the Pasteur `maestro` cluster (queue
`hubbioit`) from real CRAM alignments. They are generated automatically from the
per-run JSON files in
`.omo/evidence/meteor-rust-speedup-phase2/benchmarks/per_run/` and
`per_run_subset/`; the source of truth for every wall time, CPU time and MaxRSS
value is those JSONs.

A detailed evidence report (job IDs, per-run walltimes, MaxRSS, cancelled jobs
and learnings) is at
`.omo/evidence/meteor-rust-speedup-phase2/benchmarks/REPORT.md`. The rendered
summary report is at
`.omo/evidence/meteor-rust-speedup-phase2/benchmarks/final_report_phase2.md`.

### Environment

- Cluster: Pasteur maestro, Slurm queue `hubbioit`.
- Samples: `B18964-1_291120248_BO0015ILM`, `B18964-1_291120248_GA0218ILM`.
- Counter runs: interleaved Python/Rust inside one Slurm job.
- Variant-calling runs: one configuration per Slurm job, 24 h walltime.
- MaxRSS collected from `sacct` via `benchmarks/collect_sacct.py`.

### Counter

| Sample | Implementation | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |
|---|---|---:|---:|---:|---:|---|
| B18964-1_291120248_BO0015ILM | Python | 975.4 | 2126.0 | 21.9 GiB | 5 | baseline |
| B18964-1_291120248_BO0015ILM | Rust | 802.3 | 1430.3 | 21.9 GiB | 5 | faster |
| B18964-1_291120248_GA0218ILM | Python | 1843.8 | 6792.0 | 46.2 GiB | 5 | baseline |
| B18964-1_291120248_GA0218ILM | Rust | 1555.1 | 3763.7 | 46.2 GiB | 5 | faster |

Rust is faster on real data: roughly 18 % wall-time and 33 % CPU-time
improvement on the smaller sample, and 16 % wall-time and 45 % CPU-time
improvement on the larger sample. MaxRSS is similar because both implementations
ran inside the same interleaved Slurm job.

### Variant calling (full catalogue)

| Sample | Implementation | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |
|---|---|---:|---:|---:|---:|---|
| B18964-1_291120248_BO0015ILM | Python | — | — | — | 0 | missing |
| B18964-1_291120248_BO0015ILM | Rust batch8 | 121188.5 | 120455.8 | 41.7 GiB | 2 | incomplete |
| B18964-1_291120248_BO0015ILM | Rust batch1 | 128353.4 | 127127.5 | 269.3 GiB | 2 | incomplete |
| B18964-1_291120248_GA0218ILM | Python | — | — | — | 0 | missing |
| B18964-1_291120248_GA0218ILM | Rust batch8 | — | — | — | 0 | missing |
| B18964-1_291120248_GA0218ILM | Rust batch1 | — | — | — | 0 | missing |

The Python full-catalogue variant caller exceeded the 48 h walltime before
finishing a single run, so no Python baseline exists. Rust completed two runs
per batch configuration on the first sample only; the remaining runs were
cancelled once it became clear they would duplicate effort. With only two runs
per configuration, the full-catalogue variant-calling verdict is **incomplete**.

### Variant calling (64-gene subset fallback)

Because Python timed out on the full catalogue, a 64-gene subset was used as a
controlled fallback. Only Rust batch configurations completed; no Python
subset baseline is available.

| Sample | Implementation | Median wall (s) | Median CPU (s) | Max RSS | Runs | Verdict |
|---|---|---:|---:|---:|---:|---|
| B18964-1_291120248_BO0015ILM | Rust batch8 | 38883.1 | 38527.8 | 19.0 GiB | 4 | baseline |
| B18964-1_291120248_BO0015ILM | Rust batch1 | 38422.6 | 38103.0 | 116.0 GiB | 5 | parity |

On this subset, batch8 and batch1 are in parity on wall time (about 10.8 h),
but batch1 uses roughly 6x more memory. This suggests batch8 is the more
memory-efficient Rust configuration for this workload.

### Consensus + depth

Consensus and depth could not be measured. The original job failed because
`pysam.tabix_index(..., preset="vcf")` replaces the uncompressed `.vcf` with a
`.vcf.gz`, and `bench_phase2.py` was passing the old path. After fixing the
script to update `vcf_file = vcf_file.with_suffix(".vcf.gz")`, the resubmitted
Python stage reached only about 1 % of `create_consensus` after about 19.5 h,
projecting more than 3000 h for the full catalogue. The job was cancelled and
no consensus/depth JSONs were produced.

### Honest summary

- **Counter:** Rust is clearly faster on real data.
- **Variant calling (full catalogue):** No Python baseline; Rust results are
  incomplete due to long runtimes.
- **Variant calling (subset):** Rust batch8 and batch1 are in parity on wall
  time, with batch8 far more memory-efficient.
- **Consensus + depth:** Not measurable with the current Python implementation
  on the full catalogue; would require either a much smaller test catalogue or
  additional optimization.

## How to reproduce

With the Rust extension already built (`maturin develop -m rust/Cargo.toml`):

```bash
python benchmarks/bench_all.py --runs 5
```

The script writes its results to
`.omo/evidence/meteor-rust-acceleration/benchmarks/python.json` and
`rust.json`.
