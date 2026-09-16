#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pandas",
#     "pysam",
#     "typer",
# ]
# ///

"""Real-data benchmark driver for phase-2 Rust acceleration.

Runs one benchmark configuration per invocation and writes a single per-run JSON
file. Designed to be driven by Slurm jobs so long-running variant-calling runs
can be parallelised, while short counter/consensus+depth runs can be interleaved
inside one process.

No real data paths are hard-coded; pass them via CLI options or environment
variables.
"""

from __future__ import annotations

import csv
import json
import logging
import lzma
import os
import statistics
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NoReturn

import pandas as pd
import pysam
import subprocess
import typer

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from meteor.counter import Counter
from meteor.session import Component
from meteor.variantcalling import VariantCalling

try:
    import meteor_core
except ImportError:
    meteor_core = None  # type: ignore[assignment]


@dataclass
class RunMeasurement:
    """One timed execution."""

    wall_seconds: float
    cpu_seconds: float
    timestamp: str
    job_id: str | None = None
    max_rss_kb: int | None = None


@dataclass
class BenchmarkConfig:
    """Configuration metadata recorded with every run."""

    component: str
    sample: str
    implementation: str
    threads: int
    env: dict[str, str | None] = field(default_factory=dict)
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class BenchmarkRun:
    """Root object written to a per-run JSON file."""

    config: BenchmarkConfig
    run_index: int
    measurement: RunMeasurement
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "config": asdict(self.config),
            "run_index": self.run_index,
            "measurement": asdict(self.measurement),
            "notes": self.notes,
        }


def _fatal(message: str) -> NoReturn:
    typer.echo(message, err=True)
    raise typer.Exit(1)


def _reference_fasta(ref_dir: Path, tmp_path: Path) -> Path:
    """Return an uncompressed, indexed reference FASTA in tmp_path."""
    import gzip

    candidate = ref_dir / "reference.fa"
    if candidate.is_symlink() or candidate.is_file():
        src = candidate
    else:
        src = next(ref_dir.rglob("*.fa*"))
    src = src.resolve()
    out_fa = tmp_path / "reference.fa"
    out_fai = Path(f"{out_fa}.fai")

    with src.open("rb") as fh:
        magic = fh.read(2)
    is_compressed = src.suffix == ".gz" or magic == b"\x1f\x8b"
    if is_compressed:
        with gzip.open(src, "rb") as reader:
            out_fa.write_bytes(reader.read())
    else:
        out_fa.write_bytes(src.read_bytes())

    if not out_fai.is_file():
        pysam.faidx(str(out_fa))
    return out_fa


def _read_bed_intervals(catalogue: Path) -> list[tuple[str, int, int]]:
    intervals: list[tuple[str, int, int]] = []
    with catalogue.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 3:
                continue
            intervals.append((row[0], int(row[1]), int(row[2])))
    return intervals


def _stage1_and_ref(
    mapping_dir: Path, ref_dir: Path
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    stage1_json = next(mapping_dir.rglob("*_census_stage_1.json"))
    stage1_data = json.loads(stage1_json.read_text(encoding="utf-8"))
    ref_json_path = next(ref_dir.glob("*_reference.json"))
    ref_json = json.loads(ref_json_path.read_text(encoding="utf-8"))
    return stage1_json, stage1_data, ref_json


def _write_run(out_dir: Path, run: BenchmarkRun) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    slug = (
        f"{run.config.component}_"
        f"{run.config.sample}_"
        f"{run.config.implementation}_"
        f"run{run.run_index:02d}.json"
    )
    path = out_dir / slug
    path.write_text(json.dumps(run.to_dict(), indent=2), encoding="utf-8")
    return path


def _measure_once(func: callable[[], Any]) -> RunMeasurement:
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    func()
    wall_seconds = time.perf_counter() - wall_start
    cpu_seconds = time.process_time() - cpu_start
    return RunMeasurement(
        wall_seconds=wall_seconds,
        cpu_seconds=cpu_seconds,
        timestamp=datetime.now(timezone.utc).isoformat(),
        job_id=os.environ.get("SLURM_JOB_ID"),
    )


def _build_component(
    mapping_dir: Path, ref_dir: Path, tmp_path: Path, threads: int
) -> type[Component]:
    meteor = Component(
        threads=threads,
        fastq_dir=mapping_dir,
        mapping_dir=mapping_dir,
        ref_dir=ref_dir,
        tmp_path=tmp_path,
    )
    return meteor


def _run_counter(
    mapping_dir: Path,
    ref_dir: Path,
    sample_name: str,
    counting_type: str,
    use_rust: bool,
    threads: int,
    run_index: int,
    out_dir: Path,
) -> Path:
    stage1_json, stage1_data, ref_json = _stage1_and_ref(mapping_dir, ref_dir)
    identity_threshold = stage1_data.get("counting", {}).get(
        "identity_threshold", ref_json.get("reference_info", {}).get("identity_threshold", 0.95)
    )
    raw_cram = mapping_dir / stage1_data["mapping"]["mapping_file"]
    if not raw_cram.is_file():
        raw_cram = mapping_dir / f"{sample_name}_raw.cram"

    with tempfile.TemporaryDirectory(prefix=f"bench_counter_{counting_type}_") as tmp_raw:
        tmp_path = Path(tmp_raw)
        meteor = _build_component(mapping_dir, ref_dir, tmp_path, threads)
        meteor.use_rust_counter = use_rust
        counter = Counter(
            meteor=meteor,
            counting_type=counting_type,
            mapping_type="end-to-end",
            trim=80,
            identity_user=None,
            alignment_number=100,
            core_size=100,
            keep_all_alignments=False,
            keep_filtered_alignments=False,
            json_data={},
            identity_threshold=identity_threshold,
        )
        stage1_copy = tmp_path / stage1_json.name
        stage1_copy.write_text(json.dumps(stage1_data), encoding="utf-8")
        count_file = tmp_path / f"{sample_name}.tsv.xz"
        filtered_cram = tmp_path / f"{sample_name}.cram"

        measurement = _measure_once(
            lambda: counter.launch_counting(
                raw_cram,
                filtered_cram,
                count_file,
                ref_json,
                json.loads(stage1_copy.read_text(encoding="utf-8")),
                stage1_copy,
            )
        )

    config = BenchmarkConfig(
        component="counter",
        sample=sample_name,
        implementation="rust" if use_rust else "python",
        threads=threads,
        env={
            "METEOR_USE_RUST_COUNTER": "1" if use_rust else "0",
        },
        extra={"counting_type": counting_type},
    )
    run = BenchmarkRun(
        config=config,
        run_index=run_index,
        measurement=measurement,
    )
    return _write_run(out_dir, run)


def _run_variantcalling(
    mapping_dir: Path,
    ref_dir: Path,
    sample_name: str,
    use_rust: bool,
    batch_size: int | None,
    threads: int,
    run_index: int,
    out_dir: Path,
) -> Path:
    if use_rust and meteor_core is None:
        _fatal("Rust variant calling requested but meteor_core is not available")

    stage1_json, stage1_data, ref_json = _stage1_and_ref(mapping_dir, ref_dir)
    strain_parent = out_dir / f"strain_{sample_name}_{run_index:02d}"
    strain_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"bench_vc_{sample_name}_", dir=strain_parent) as tmp_raw:
        tmp_path = Path(tmp_raw)
        strain_dir = tmp_path / "strain" / sample_name
        strain_dir.mkdir(parents=True)

        meteor = _build_component(mapping_dir, ref_dir, tmp_path, threads)
        meteor.use_rust_variant_calling = use_rust
        meteor.strain_dir = strain_dir

        census = {
            "mapped_sample_dir": mapping_dir,
            "directory": strain_dir,
            "census": stage1_data,
            "reference": ref_json,
            "Stage3FileName": strain_dir / f"{sample_name}_census_stage_3.json",
        }

        variant_caller = VariantCalling(
            meteor=meteor,
            census=census,
            max_depth=100,
            min_depth=3,
            min_snp_depth=3,
            min_frequency=0.1,
            ploidy=1,
            core_size=10,
        )

        env_before = os.environ.get("METEOR_FREEBAYES_BATCH_SIZE")
        if batch_size is not None:
            os.environ["METEOR_FREEBAYES_BATCH_SIZE"] = str(batch_size)

        try:
            measurement = _measure_once(variant_caller.execute)
        finally:
            if env_before is None:
                os.environ.pop("METEOR_FREEBAYES_BATCH_SIZE", None)
            else:
                os.environ["METEOR_FREEBAYES_BATCH_SIZE"] = env_before

    impl = "python"
    if use_rust:
        impl = f"rust-batch{batch_size}" if batch_size else "rust-default"

    config = BenchmarkConfig(
        component="variantcalling",
        sample=sample_name,
        implementation=impl,
        threads=threads,
        env={
            "METEOR_USE_RUST_VARIANT_CALLING": "1" if use_rust else "0",
            "METEOR_FREEBAYES_BATCH_SIZE": str(batch_size) if batch_size is not None else "default",
        },
    )
    run = BenchmarkRun(
        config=config,
        run_index=run_index,
        measurement=measurement,
    )
    return _write_run(out_dir, run)


def _make_subset_ref_dir(
    ref_dir: Path, catalogue: Path | None, tmp_path: Path, n_genes: int
) -> Path:
    """Create a temporary reference directory containing only the first N genes.

    The full MSP map and catalogue are truncated; the remaining files are
    symlinked from the full reference so that VariantCalling.execute can run on a
    small, well-defined gene set.
    """
    full_catalogue = catalogue if catalogue is not None else ref_dir / "catalogue.bed"
    if not full_catalogue.is_file():
        _fatal(f"Cannot find catalogue.bed at {full_catalogue}")

    subset_dir = tmp_path / "subset_ref"
    subset_dir.mkdir(parents=True)
    db_dir = subset_dir / "database"
    db_dir.mkdir()

    with full_catalogue.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = [row for row in reader if len(row) >= 3]
    subset_rows = rows[:n_genes]
    gene_ids = {row[0] for row in subset_rows}

    subset_catalogue = subset_dir / "catalogue.bed"
    with subset_catalogue.open("w") as fh:
        for row in subset_rows:
            fh.write("\t".join(row[:3]) + "\n")

    msp_map = ref_dir / "database" / "msp_map.tsv"
    if not msp_map.is_file():
        msp_map = ref_dir / "msp_map.tsv"
    if msp_map.is_file():
        with msp_map.open() as fh:
            header = fh.readline()
            cols = header.strip().split("\t")
            try:
                gene_id_idx = cols.index("gene_id")
            except ValueError:
                gene_id_idx = 1
            msp_rows = [header] + [
                line for line in fh if len(line.split("\t")) > gene_id_idx and line.split("\t")[gene_id_idx].strip() in gene_ids
            ]
        (db_dir / "msp_map.tsv").write_text("".join(msp_rows), encoding="utf-8")
    else:
        _fatal(f"Cannot find msp_map.tsv in {ref_dir}/database")

    ref_json_path = next(ref_dir.glob("*_reference.json"))
    ref_json = json.loads(ref_json_path.read_text(encoding="utf-8"))
    ref_json["reference_file"]["database_dir"] = "database"
    ref_json["reference_file"]["fasta_dir"] = "."
    (subset_dir / ref_json_path.name).write_text(json.dumps(ref_json), encoding="utf-8")

    # Build a small FASTA containing only the selected genes so downstream tools
    # do not load the whole catalogue.  The reference is block-gzip compressed,
    # so re-compress the subset and index it for pysam/freebayes.
    for fasta in ref_dir.glob("*.fa*"):
        if fasta.suffix in {".fai", ".gzi"}:
            continue
        subset_fasta = subset_dir / fasta.name
        tmp_plain = subset_fasta.with_suffix(subset_fasta.suffix + ".plain")
        with pysam.FastaFile(str(fasta)) as fh_in, tmp_plain.open("w") as fh_out:
            for ref_name in fh_in.references:
                if ref_name in gene_ids:
                    seq = fh_in.fetch(ref_name)
                    fh_out.write(f">{ref_name}\n{seq}\n")
        subprocess.run(
            ["bgzip", "-f", "-c", str(tmp_plain)],
            check=True,
            stdout=subset_fasta.open("wb"),
        )
        tmp_plain.unlink(missing_ok=True)
        pysam.faidx(str(subset_fasta))

    for feather in (ref_dir / "database").glob("*.feather"):
        if "msp" not in feather.name.lower():
            (db_dir / feather.name).symlink_to(feather.resolve())

    return subset_dir


def _run_variantcalling_subset(
    mapping_dir: Path,
    ref_dir: Path,
    catalogue: Path | None,
    sample_name: str,
    use_rust: bool,
    batch_size: int | None,
    n_genes: int,
    threads: int,
    run_index: int,
    out_dir: Path,
) -> Path:
    if use_rust and meteor_core is None:
        _fatal("Rust variant calling requested but meteor_core is not available")

    stage1_json, stage1_data, _ref_json = _stage1_and_ref(mapping_dir, ref_dir)
    strain_parent = out_dir / f"strain_{sample_name}_{run_index:02d}"
    strain_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"bench_vc_subset_{sample_name}_", dir=strain_parent) as tmp_raw:
        tmp_path = Path(tmp_raw)
        subset_ref = _make_subset_ref_dir(ref_dir, catalogue, tmp_path, n_genes)
        strain_dir = tmp_path / "strain" / sample_name
        strain_dir.mkdir(parents=True)

        meteor = _build_component(mapping_dir, subset_ref, tmp_path, threads)
        meteor.use_rust_variant_calling = use_rust
        meteor.strain_dir = strain_dir

        ref_json = json.loads(next(subset_ref.glob("*_reference.json")).read_text(encoding="utf-8"))
        census = {
            "mapped_sample_dir": mapping_dir,
            "directory": strain_dir,
            "census": stage1_data,
            "reference": ref_json,
            "Stage3FileName": strain_dir / f"{sample_name}_census_stage_3.json",
        }

        variant_caller = VariantCalling(
            meteor=meteor,
            census=census,
            max_depth=100,
            min_depth=3,
            min_snp_depth=3,
            min_frequency=0.1,
            ploidy=1,
            core_size=10,
        )

        env_before = os.environ.get("METEOR_FREEBAYES_BATCH_SIZE")
        if batch_size is not None:
            os.environ["METEOR_FREEBAYES_BATCH_SIZE"] = str(batch_size)

        try:
            measurement = _measure_once(variant_caller.execute)
        finally:
            if env_before is None:
                os.environ.pop("METEOR_FREEBAYES_BATCH_SIZE", None)
            else:
                os.environ["METEOR_FREEBAYES_BATCH_SIZE"] = env_before

    impl = "python"
    if use_rust:
        impl = f"rust-batch{batch_size}" if batch_size else "rust-default"

    config = BenchmarkConfig(
        component="variantcalling-subset",
        sample=sample_name,
        implementation=impl,
        threads=threads,
        env={
            "METEOR_USE_RUST_VARIANT_CALLING": "1" if use_rust else "0",
            "METEOR_FREEBAYES_BATCH_SIZE": str(batch_size) if batch_size is not None else "default",
        },
        extra={"subset_genes": n_genes},
    )
    run = BenchmarkRun(
        config=config,
        run_index=run_index,
        measurement=measurement,
        notes=[f"limited to first {n_genes} genes from catalogue"],
    )
    return _write_run(out_dir, run)


def _run_consensus_depth(
    mapping_dir: Path,
    ref_dir: Path,
    catalogue: Path,
    sample_name: str,
    use_rust: bool,
    threads: int,
    run_index: int,
    out_dir: Path,
) -> Path:
    if use_rust and meteor_core is None:
        _fatal("Rust consensus+depth requested but meteor_core is not available")

    _, stage1_data, ref_json = _stage1_and_ref(mapping_dir, ref_dir)
    cram = mapping_dir / f"{sample_name}.cram"
    if not cram.is_file():
        cram = mapping_dir / stage1_data["mapping"]["mapping_file"]

    with tempfile.TemporaryDirectory(prefix=f"bench_cd_{sample_name}_") as tmp_raw:
        tmp_path = Path(tmp_raw)
        ref_fa = _reference_fasta(ref_dir, tmp_path)
        intervals = _read_bed_intervals(catalogue)
        bed = tmp_path / "regions.bed"
        with bed.open("w") as fh:
            for gene_id, start, end in intervals:
                fh.write(f"{gene_id}\t{start}\t{end}\n")

        matrix_rows = [
            {"gene_id": gene_id, "gene_length": end - start, "coverage": 1000}
            for gene_id, start, end in intervals
        ]
        matrix_file = tmp_path / "matrix.tsv.xz"
        pd.DataFrame(matrix_rows).to_csv(
            matrix_file, sep="\t", index=False, header=True, compression="xz"
        )

        # Generate a VCF once per run. The VCF generation itself is not timed;
        # the benchmark measures only the consensus+depth hot paths that depend on
        # the Python vs Rust implementations.
        vcf_file = tmp_path / "consensus.vcf"
        if meteor_core is not None:
            options = meteor_core.FreebayesOptions(
                min_snp_depth=1,
                min_frequency=0.1,
                ploidy=1,
            )
            meteor_core.call_variants_parallel(
                str(cram.resolve()),
                str(ref_fa.resolve()),
                str(bed),
                "freebayes",
                options,
                n_threads=threads,
                output_path=str(vcf_file),
            )
            pysam.tabix_index(str(vcf_file), preset="vcf")
            # pysam.tabix_index compresses the VCF and removes the uncompressed
            # file; point to the generated .vcf.gz so create_consensus can open it.
            vcf_file = vcf_file.with_suffix(".vcf.gz")
        else:
            _fatal("meteor_core is required for consensus_depth benchmark VCF generation")

        meteor = _build_component(mapping_dir, ref_dir, tmp_path, threads)
        meteor.use_rust_variant_calling = use_rust
        vc = VariantCalling(
            meteor=meteor,
            census={},
            max_depth=8000,
            min_depth=3,
            min_snp_depth=1,
            min_frequency=0.1,
            ploidy=1,
            core_size=100,
        )
        vc.matrix_file = matrix_file

        def _work() -> None:
            low_cov_sites, gene_ignore = vc.filter_low_cov_sites(cram, ref_fa)
            consensus_file = tmp_path / "consensus.fasta.xz"
            vc.create_consensus(ref_fa, consensus_file, low_cov_sites, gene_ignore, vcf_file, bed)

        measurement = _measure_once(_work)

    config = BenchmarkConfig(
        component="consensus_depth",
        sample=sample_name,
        implementation="rust" if use_rust else "python",
        threads=threads,
        env={
            "METEOR_USE_RUST_VARIANT_CALLING": "1" if use_rust else "0",
        },
    )
    run = BenchmarkRun(
        config=config,
        run_index=run_index,
        measurement=measurement,
    )
    return _write_run(out_dir, run)


app = typer.Typer(add_completion=False, pretty_exceptions_short=True)


@app.command()
def counter(
    mapping_dir: Path = typer.Argument(..., help="Directory containing mapped sample CRAM/JSON."),
    ref_dir: Path = typer.Argument(..., help="Directory containing reference JSON."),
    sample_name: str = typer.Argument(..., help="Sample name used in output filenames."),
    out_dir: Path = typer.Option(..., "--out-dir", help="Directory for per-run JSON files."),
    counting_type: str = typer.Option("smart_shared", "--counting-type"),
    use_rust: bool = typer.Option(False, "--use-rust/--python"),
    threads: int = typer.Option(16, "--threads"),
    run_index: int = typer.Option(0, "--run-index"),
) -> None:
    """Benchmark one counter configuration."""
    path = _run_counter(
        mapping_dir, ref_dir, sample_name, counting_type, use_rust, threads, run_index, out_dir
    )
    typer.echo(path)


@app.command()
def variantcalling(
    mapping_dir: Path = typer.Argument(..., help="Directory containing mapped sample CRAM/JSON."),
    ref_dir: Path = typer.Argument(..., help="Directory containing reference JSON."),
    sample_name: str = typer.Argument(..., help="Sample name used in output filenames."),
    out_dir: Path = typer.Option(..., "--out-dir", help="Directory for per-run JSON files."),
    use_rust: bool = typer.Option(False, "--use-rust/--python"),
    batch_size: int | None = typer.Option(None, "--batch-size"),
    threads: int = typer.Option(16, "--threads"),
    run_index: int = typer.Option(0, "--run-index"),
) -> None:
    """Benchmark one variant-calling configuration."""
    path = _run_variantcalling(
        mapping_dir, ref_dir, sample_name, use_rust, batch_size, threads, run_index, out_dir
    )
    typer.echo(path)


@app.command()
def variantcalling_subset(
    mapping_dir: Path = typer.Argument(..., help="Directory containing mapped sample CRAM/JSON."),
    ref_dir: Path = typer.Argument(..., help="Directory containing reference JSON."),
    sample_name: str = typer.Argument(..., help="Sample name used in output filenames."),
    out_dir: Path = typer.Option(..., "--out-dir", help="Directory for per-run JSON files."),
    catalogue: Path | None = typer.Option(None, "--catalogue", help="Catalogue BED file."),
    use_rust: bool = typer.Option(False, "--use-rust/--python"),
    batch_size: int | None = typer.Option(None, "--batch-size"),
    n_genes: int = typer.Option(64, "--n-genes"),
    threads: int = typer.Option(8, "--threads"),
    run_index: int = typer.Option(0, "--run-index"),
) -> None:
    """Benchmark variant calling on the first N genes of the catalogue."""
    path = _run_variantcalling_subset(
        mapping_dir, ref_dir, catalogue, sample_name, use_rust, batch_size, n_genes, threads, run_index, out_dir
    )
    typer.echo(path)


@app.command()
def consensus_depth(
    mapping_dir: Path = typer.Argument(..., help="Directory containing mapped sample CRAM/JSON."),
    ref_dir: Path = typer.Argument(..., help="Directory containing reference JSON."),
    catalogue: Path = typer.Argument(..., help="Catalogue BED file with gene intervals."),
    sample_name: str = typer.Argument(..., help="Sample name used in output filenames."),
    out_dir: Path = typer.Option(..., "--out-dir", help="Directory for per-run JSON files."),
    use_rust: bool = typer.Option(False, "--use-rust/--python"),
    threads: int = typer.Option(16, "--threads"),
    run_index: int = typer.Option(0, "--run-index"),
) -> None:
    """Benchmark one consensus+depth configuration."""
    path = _run_consensus_depth(
        mapping_dir, ref_dir, catalogue, sample_name, use_rust, threads, run_index, out_dir
    )
    typer.echo(path)


@app.command()
def interleaved(
    mapping_dir: Path = typer.Argument(...),
    ref_dir: Path = typer.Argument(...),
    sample_name: str = typer.Argument(...),
    out_dir: Path = typer.Option(..., "--out-dir"),
    component: str = typer.Option("counter", "--component"),
    runs: int = typer.Option(5, "--runs"),
    start_run: int = typer.Option(0, "--start-run"),
    threads: int = typer.Option(16, "--threads"),
    catalogue: Path | None = typer.Option(None, "--catalogue"),
) -> None:
    """Run interleaved Python/Rust configs for the given component.

    For the counter component, this alternates Python/Rust for both
    smart_shared and total counting types. For consensus_depth, it alternates
    Python and Rust implementations.
    """
    if component == "counter":
        configs: list[tuple[str, bool, dict[str, Any]]] = [
            ("smart_shared", False, {"counting_type": "smart_shared"}),
            ("smart_shared", True, {"counting_type": "smart_shared"}),
            ("total", False, {"counting_type": "total"}),
            ("total", True, {"counting_type": "total"}),
        ]
        for run_index in range(start_run, start_run + runs):
            for impl, use_rust, extra in configs:
                _run_counter(
                    mapping_dir,
                    ref_dir,
                    sample_name,
                    extra["counting_type"],
                    use_rust,
                    threads,
                    run_index,
                    out_dir,
                )
                typer.echo(
                    f"counter {impl} run {run_index} done: "
                    f"{sample_name} {'rust' if use_rust else 'python'}"
                )
    elif component == "consensus_depth":
        if catalogue is None:
            _fatal("--catalogue is required for consensus_depth interleaved runs")
        for run_index in range(start_run, start_run + runs):
            for use_rust in (False, True):
                _run_consensus_depth(
                    mapping_dir,
                    ref_dir,
                    catalogue,
                    sample_name,
                    use_rust,
                    threads,
                    run_index,
                    out_dir,
                )
                typer.echo(
                    f"consensus_depth run {run_index} done: "
                    f"{'rust' if use_rust else 'python'}"
                )
    else:
        _fatal(f"Unknown component for interleaved mode: {component}")


if __name__ == "__main__":
    app()
