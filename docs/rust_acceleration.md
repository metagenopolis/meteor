# Rust acceleration

Meteor ships with an **optional** Rust extension (`meteor-core`) that accelerates
CRAM counting and strain-profiling hot paths. The extension is built with
[maturin](https://www.maturin.rs/) and is **never required** at runtime: every
Rust path has a pure-Python fallback.

## What is accelerated

| Meteor step | Rust helper | Python fallback |
|-------------|-------------|-----------------|
| `meteor mapping` / `meteor profile` counting | `meteor_core.count_msp` | `meteor.counter.Counter.launch_counting` |
| `meteor strain` read pileup | `meteor_core.count_reads_in_gene` | `meteor.variantcalling.VariantCalling.count_reads_in_gene` |
| `meteor strain` consensus building | `meteor_core.create_consensus` | `meteor.variantcalling.VariantCalling.create_consensus` |
| `meteor strain` freebayes dispatch | `meteor_core.call_variants_parallel` | `concurrent.futures.ProcessPoolExecutor` |
| Streaming CRAM records | `meteor_core.stream_cram_records` | `meteor_core.cram_records` (buffered) |
| Output helpers | `meteor_core.write_vcf_text`, `bgzip_file`, `write_bcf` | Python gzip / pysam equivalents |

## Build the extension

From the repository root:

```bash
pip install maturin
maturin develop -m rust/Cargo.toml
```

For a release wheel:

```bash
maturin build -m rust/Cargo.toml --release
```

The Python package continues to use `poetry-core` as its default build backend.
To install with the optional Rust extra (once the `meteor-core` wheel is
published):

```bash
pip install meteor[rust]
```

## Enable the extension

You can enable Rust acceleration with CLI flags or environment variables.

### Counting

```bash
meteor mapping ... --use-rust-counter
METEOR_USE_RUST_COUNTER=1 meteor mapping ...
```

### Variant calling

```bash
meteor strain ... --use-rust-variant-calling
METEOR_USE_RUST_VARIANT_CALLING=1 meteor strain ...
```

The old environment variable `METEOR_USE_RUST_VARIANT=1` is still accepted but
logs a deprecation warning.

### Combined flag

```bash
# Enables every Rust helper supported by the current command
meteor mapping ... --use-rust
meteor strain ... --use-rust
METEOR_USE_RUST=1 meteor mapping ...
METEOR_USE_RUST=1 meteor strain ...
```

`--use-rust` / `METEOR_USE_RUST=1` is equivalent to setting both the counting
and variant-calling flags for the command it is applied to.

## Runtime fallback

If `meteor-core` is not installed **or** a Rust helper raises any runtime
error (including a Rust panic), Meteor logs a warning and silently falls back
to the Python implementation. The pipeline is never aborted just because the
optional extension is missing or failed.

You can verify which implementation ran by checking the log output:

```
INFO: Used Rust counter implementation
WARNING: Rust counter requested but meteor_core is not available. Falling back to Python.
```

## Continuous integration

The GitHub Actions workflow builds the Rust extension on Ubuntu and macOS and
runs `cargo clippy -- -D warnings`. The Python test matrix also builds the
extension before running pytest so the Rust parity tests are exercised on
CI.

## Performance caveats

See [docs/benchmarks.md](benchmarks.md) for the measured numbers. The report is
honest about both wins and regressions:

- **Counter**: the Rust `meteor_core.count_msp` implementation is currently on
  par with or slightly slower than the Python hot loop on the repository
  fixture. The reason is that `count_msp` returns per-gene **read lists**
  (Python string objects for every classified read) across the PyO3 boundary,
  so the conversion cost dominates the small-fixture run time. A future redesign
  that returns only per-gene aggregates from Rust is the expected path to a
  real speed-up.
- **Variant calling**: the Rust-accelerated `meteor strain` path is slower on
  the small fixture because most of the time is spent launching short
  `freebayes` subprocesses. The Rust helpers are expected to become beneficial
  on larger catalogues with many more regions, but that has not been
  benchmarked yet.

## Benchmarks

See [docs/benchmarks.md](benchmarks.md) for replicated benchmark results on the
repository fixtures. The report includes both speedups and cases where the
Rust path is slower on small inputs, so expectations are honest.

## Conda note

The Bioconda `meteor` package is `noarch: python` and does not ship the
compiled Rust extension. Conda users who want Rust acceleration should install
or build `meteor-core` separately.
