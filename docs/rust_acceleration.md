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
  so the conversion cost dominates the small-fixture run time. Phase 2
  implements exactly the redesign described here, see
  [Phase 2 extensions](#phase-2-extensions).
- **Variant calling**: the Rust-accelerated `meteor strain` path is slower on
  the small fixture because most of the time is spent launching short
  `freebayes` subprocesses. Phase 2 batches those invocations, see
  [Phase 2 extensions](#phase-2-extensions).

## Benchmarks

See [docs/benchmarks.md](benchmarks.md) for replicated benchmark results on the
repository fixtures. The report includes both speedups and cases where the
Rust path is slower on small inputs, so expectations are honest.

## Conda note

The Bioconda `meteor` package is `noarch: python` and does not ship the
compiled Rust extension. Conda users who want Rust acceleration should install
or build `meteor-core` separately.

## Phase 2 extensions

Phase 2 (branch `meteor-rust-speedup-phase2`, stacked on PR #126) removes the
overheads measured in phase 1. The counting hot path no longer ships per-read
data across the Rust/Python boundary, the whole count table can be written by
Rust, freebayes is launched far fewer times, and the per-gene coverage pileup
becomes a single pass over the CRAM. Everything stays opt-in: the same
`--use-rust*` flags as before select the Rust paths, and every Rust failure
still falls back to Python with a logged warning.

Measured results live in [docs/benchmarks.md](benchmarks.md). This section
documents the new APIs, the environment knobs and the go/no-go decisions that
shipped them.

### Aggregates-only counter APIs

`meteor_core` gains two counting entry points:

| Function | Arguments | Returns |
|----------|-----------|---------|
| `count_msp_aggregates` | `(cram_path, msp_map_path, identity_threshold, counting_type)` | list of `AggregateRow` |
| `count_msp_write_tsv` | `(cram_path, msp_map_path, out_tsv_path, identity_threshold, counting_type)` | `(rows_written, counted_reads)` |

- `AggregateRow` carries `msp`, `gene`, `count` and `reads` for one gene.
  `reads` is a single newline-delimited string holding the sorted,
  deduplicated read ids that contributed to the gene, so the Rust/Python
  boundary is crossed roughly once per gene instead of once per read. The
  Python side derives `counted_reads` from the distinct read names across
  rows, converts whole-number counts back to ints
  (`Counter._normalise_count_value`), rebuilds the gene-length table from the
  CRAM header and writes the TSV with `Counter.write_stat`.
- `count_msp_write_tsv` writes the complete count table from Rust: xz
  compression preset 0, the same `gene_id`, `gene_length`, `value` column
  header, rows sorted by gene id, and CPython-compatible float formatting. The
  output is byte-identical to `Counter.write_stat`, and only the output path
  crosses the boundary (no per-gene data at all). Python only updates the
  stage-1 JSON report afterwards.
- Both functions take the **MSP map file path** instead of the reference
  FASTA. The reference is located next to it as `reference.fa`,
  `reference.fa.gz`, `reference.fasta` or `reference.fasta.gz`.
- `counting_type` accepts `smart_shared`, `unique` and `total`, with the same
  semantics as `count_msp`. The phase-1 `count_msp` API is unchanged and
  remains the parity reference.

The dispatcher `Counter._launch_counting_rust` walks a capability chain, so an
extension build that predates one of the new functions still works:

1. `count_msp_write_tsv` when present,
2. `count_msp_aggregates` when present,
3. `count_msp` (phase-1 API),
4. otherwise a `RuntimeError`, which the caller turns into the usual warning
   and Python fallback.

Each step checks for the function with `hasattr` and catches any exception:
a failure logs a warning and drops to the next function. The Rust counter is
only attempted when `--use-rust-counter` (or `METEOR_USE_RUST_COUNTER=1` /
`METEOR_USE_RUST=1`) is set and `--kf` is not, because keeping filtered
alignments requires the per-read data only the Python path can emit.

### Single-pass coverage pileup

`meteor_core.depth_per_gene(cram_path, ref_path, genes, max_depth)` replaces
the per-gene `fetch` + `pileup` loop for the low-coverage detection step.
`genes` is a list of `(gene_id, start, end)` 0-based half-open intervals and
the call returns one `GeneDepth { gene, depths }` per input interval, in input
order.

- One streaming pass over the CRAM accumulates every gene at once. Intervals
  are grouped per contig and sorted by start, and each coverage block is
  applied with a binary search plus a sweep pointer, so the cost scales with
  the records rather than records times genes.
- Semantics match the phase-1 pileup: positions covered by a deletion or a
  reference skip are not counted, depths are capped at `max_depth`, positions
  with zero coverage are present as `0`, and overlapping gene intervals
  accumulate independently.
- Unmapped, secondary, duplicate and QC-failed records are skipped.
- Invalid intervals (`start >= end`) and gene ids missing from the CRAM header
  raise a `ValueError`. The CRAM must be indexed.

`VariantCalling.filter_low_cov_sites` uses it when
`--use-rust-variant-calling` (or `METEOR_USE_RUST_VARIANT_CALLING=1` /
`METEOR_USE_RUST=1`) is active, with the fallback chain preserved:

1. one `meteor_core.depth_per_gene` call for all genes whose total coverage
   reaches `min_depth`,
2. per-gene `meteor_core.count_reads_in_gene` (phase-1 API) if that fails,
3. the Python `pysam` pileup path.

### Batched freebayes execution

`meteor_core.call_variants_parallel` now groups the per-thread BED chunks into
batches before launching freebayes. Chunk boundaries are computed exactly as
before (a balanced split of the BED lines), then adjacent chunks are
concatenated into batch BED files, so the set of `-t` intervals handed to
freebayes is identical to the unbatched dispatcher. One freebayes process runs
per batch instead of per chunk, batches are distributed over the same
`std::thread::scope` worker pool, and each batch keeps the existing 3600 s
timeout and stderr capture. A failing or hung batch surfaces as
`meteor_core.FreebayesError` with the batch index, exit code and stderr.
Batch outputs are merged in batch order under a single shared header; the
Python side then sorts with `bcftools` and indexes as before.

`METEOR_FREEBAYES_BATCH_SIZE=1` reproduces the previous one-process-per-chunk
behaviour exactly.

### Threaded CRAM decoding

Every Rust CRAM open now applies an htslib thread pool through
`set_threads(n)`. The counter and pileup paths share the `open_cram` helper,
the aggregates paths use their own `open_cram_with_msp_map`, and
`count_reads_in_gene` sets the pool inline, so all of them honour the same
resolver (`rust/src/threads.rs`) driven by `METEOR_RUST_THREADS`:

| `METEOR_RUST_THREADS` | Effect |
|-----------------------|--------|
| unset | `min(8, available_parallelism())`, or 1 if the system parallelism cannot be determined |
| `0` | the default above |
| positive integer | that value, capped at `available_parallelism()` |
| non-integer (e.g. `abc`) | warning on stderr, then the default |

Only the Rust path is affected. The Python/pysam side keeps its own `threads`
setting.

### Environment variables

| Variable | Used by | Semantics |
|----------|---------|-----------|
| `METEOR_RUST_THREADS` | all Rust CRAM opens | htslib decode thread pool size, see the table above |
| `METEOR_FREEBAYES_BATCH_SIZE` | `call_variants_parallel` | chunks per freebayes process; unset means `8`; `0` means `1` with a stderr warning; a positive integer is used as-is; a non-integer logs a warning and falls back to `8` |
| `METEOR_REAL_BENCH_DATA` | `tests/test_realdata_parity.py` | must be exactly `1` to run the real-data parity suite; any other value skips the module so CI stays green |
| `METEOR_BENCH_CRAM`, `METEOR_BENCH_REF`, `METEOR_BENCH_CATALOGUE`, `METEOR_BENCH_MSP_MAP` | real-data parity tests and benchmarks | absolute paths to a filtered sample CRAM, the reference directory, the gene catalogue BED and the MSP map TSV; the values are machine-local and never committed |
| `METEOR_PARITY_VCF_GENES` (default `64`), `METEOR_PARITY_CONSENSUS_GENES` (default `64`), `METEOR_PARITY_DEPTH_GENES` (default `100`) | `tests/test_realdata_parity.py` | bound how many catalogue intervals each parity assertion covers |

`METEOR_USE_RUST_COUNTER`, `METEOR_USE_RUST_VARIANT_CALLING` and
`METEOR_USE_RUST` predate phase 2 and are described under
[Enable the extension](#enable-the-extension).

### Go/no-go decisions

Each phase-2 optimisation landed only after real-data profiling justified it,
and the real-data parity suite proved the outputs unchanged
(`.omo/evidence/meteor-rust-speedup-phase2/profiling/report_phase2.md` and
`task-8-parity/REPORT.md`).

| Component | Decision | Basis |
|-----------|----------|-------|
| Counter, TSV written in Rust | GO | Profiling on real samples showed that even the aggregates API was dominated by the Rust/Python boundary: per-gene read-name strings were materialised in Rust and re-split in Python, so writing the TSV inside Rust was the right fix. |
| Counter regression and fix | initially NO-GO, then GO | On the 10.4 million-gene `hs_10_4_gut` reference the first TSV-in-Rust implementation never finished: `HeaderView::target_names()` reallocated every target name for each record, an O(n²) blow-up over the header, compounded by a redundant header clone and per-gene read-name vectors the TSV path never used. Commit `450b411` caches the target names once and skips read-name collection on the TSV path. After the fix the Rust counter produces byte-identical TSVs on real data in both `smart_shared` and `total` modes. |
| Variant calling, freebayes batching | GO | Profiling showed the Python dispatcher spends most of its wall time waiting on short freebayes processes, so reducing the process count was justified. |
| Variant calling, single-pass pileup | GO | The same profiling showed the Rust path was CPU-saturated inside the per-gene fetch and pileup work, so replacing it with the single-pass `depth_per_gene` was justified. |
| Consensus | GO, unchanged in phase 2 | `create_consensus` stayed as phase 1 shipped it; real-data parity confirmed byte-identical FASTA output, so the Rust consensus path remains safe to use. |

### Running the real-data parity suite

`tests/test_realdata_parity.py` asserts counter TSV byte identity, VCF record
equality after `bcftools norm -f ref | bcftools sort`, consensus FASTA
identity and depth-array equality on real samples. Without the env gate the
whole module reports as skipped:

```bash
# gate off: every test in the module is skipped, CI stays green
pytest tests/ -q

# gate on: absolute, machine-local paths
METEOR_REAL_BENCH_DATA=1 \
METEOR_BENCH_CRAM=/abs/path/sample.cram \
METEOR_BENCH_REF=/abs/path/ref_for_rust \
METEOR_BENCH_CATALOGUE=/abs/path/catalogue.bed \
METEOR_BENCH_MSP_MAP=/abs/path/msp_map.tsv \
pytest tests/test_realdata_parity.py -v
```

With the gate on but a path variable missing, the affected tests skip and name
the missing variable. On a mismatch the failing test writes a unified diff
under `.omo/evidence/meteor-rust-speedup-phase2/realdata/` and fails with the
path to that diff.
