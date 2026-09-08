#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pysam",
#     "numpy",
#     "pandas",
#     "typer",
#     "bgzip",
# ]
# ///

# ─── How to run ───
# 1. Install uv (if not installed):
#      curl -LsSf https://astral.sh/uv/install.sh | sh
# 2. Run directly (no venv, no pip install needed):
#      uv run tests/data/fixtures/make_fixtures.py
# 3. With custom output directory or seed:
#      uv run tests/data/fixtures/make_fixtures.py --output-dir ./fixtures --seed 123
# ──────────────────

"""Deterministic fixture generator for the meteor Rust-acceleration benchmarks."""

from __future__ import annotations

import json
import lzma
import shutil
import subprocess
import textwrap
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import bgzip
import numpy as np
import pandas as pd
import pysam
import typer
from numpy.random import Generator, default_rng
from pysam import AlignmentFile, AlignmentHeader

if TYPE_CHECKING:
    pass

app = typer.Typer(add_completion=False, pretty_exceptions_short=True)

NUCLEOTIDES = np.array(["A", "C", "G", "T"], dtype=str)
DEFAULT_GENES = 120
DEFAULT_READS = 22_000
DEFAULT_READ_LEN = 100
DEFAULT_MSP_COUNT = 10


def _write_fasta(path: Path, sequences: dict[int, str]) -> None:
    """Write a FASTA file with 80-character wrapped sequence lines."""
    with path.open("wt", encoding="utf-8") as fh:
        for gene_id, seq in sequences.items():
            fh.write(f">{gene_id}\n")
            fh.write(textwrap.fill(seq, width=80))
            fh.write("\n")


def _make_reference(
    output_dir: Path,
    rng: Generator,
    n_genes: int,
    min_len: int = 1_000,
    max_len: int = 3_000,
) -> dict[int, str]:
    """Generate reference.fa, reference.fa.fai, reference.fa.gz and indexes."""
    ref_fa = output_dir / "reference.fa"
    ref_gz = output_dir / "reference.fa.gz"

    sequences: dict[int, str] = {}
    for gene_id in range(1, n_genes + 1):
        length = int(rng.integers(min_len, max_len + 1))
        sequences[gene_id] = "".join(rng.choice(NUCLEOTIDES, size=length))

    _write_fasta(ref_fa, sequences)
    pysam.faidx(str(ref_fa))

    # Bgzip the reference so variant calling can consume it the way meteor expects.
    with ref_fa.open("rb") as plain, ref_gz.open("wb") as compressed:
        with bgzip.BGZipWriter(compressed) as writer:
            writer.write(plain.read())
    pysam.faidx(str(ref_gz))

    return sequences


def _mutate_read(read: str, rng: Generator, sub_rate: float = 0.01, indel_rate: float = 0.005) -> str:
    """Introduce substitutions and short indels into a read."""
    chars = list(read)
    i = 0
    result: list[str] = []
    while i < len(chars):
        if rng.random() < indel_rate:
            # Short deletion: skip one base.
            if rng.random() < 0.5 and i + 1 < len(chars):
                i += 1
            else:
                # Short insertion: duplicate the current base with a random partner.
                result.append(chars[i])
                result.append(rng.choice(NUCLEOTIDES))
        elif rng.random() < sub_rate:
            result.append(rng.choice(NUCLEOTIDES))
        else:
            result.append(chars[i])
        i += 1
    return "".join(result)


def _make_reads(
    output_dir: Path,
    sequences: dict[int, str],
    rng: Generator,
    n_reads: int,
    read_len: int,
) -> dict[int, int]:
    """Generate a FASTQ file and return per-gene read counts."""
    fq_file = output_dir / "sample.fq"
    gene_ids = np.array(list(sequences.keys()), dtype=int)
    gene_lengths = np.array([len(sequences[g]) for g in gene_ids], dtype=int)
    probs = gene_lengths / gene_lengths.sum()

    gene_counts: dict[int, int] = {gene_id: 0 for gene_id in gene_ids}

    with fq_file.open("wt", encoding="utf-8") as fh:
        for read_idx in range(1, n_reads + 1):
            gene_id = int(rng.choice(gene_ids, p=probs))
            seq = sequences[gene_id]
            max_start = max(1, len(seq) - read_len + 1)
            start = int(rng.integers(0, max_start))
            fragment = seq[start : start + read_len]
            mutated = _mutate_read(fragment, rng)
            quality = "".join(chr(int(rng.integers(30, 42)) + 33) for _ in mutated)
            fh.write(f"@{read_idx}\n{mutated}\n+\n{quality}\n")
            gene_counts[gene_id] += 1

    return gene_counts


def _sanitize_cram_header(cram_file: Path, ref_fa: Path) -> None:
    """Strip absolute paths and @PG records so CRAMs are byte-identical across runs."""
    tmp_cram = cram_file.with_suffix(".tmp.cram")
    with AlignmentFile(str(cram_file), "rc", reference_filename=str(ref_fa)) as infile:
        header_dict = infile.header.to_dict()
        header_dict.pop("PG", None)
        for sq in header_dict.get("SQ", []):
            sq.pop("UR", None)
        clean_header = AlignmentHeader.from_dict(header_dict)
        with AlignmentFile(
            str(tmp_cram),
            "wc",
            header=clean_header,
            reference_filename=str(ref_fa),
        ) as outfile:
            for read in infile:
                outfile.write(read)
    tmp_cram.replace(cram_file)


def _align_reads(
    output_dir: Path,
    ref_fa: Path,
    fq_file: Path,
    seed: int,
) -> Path:
    """Align reads with bowtie2 and produce a coordinate-sorted CRAM + index."""
    index_dir = output_dir / "bt2_index"
    index_dir.mkdir(exist_ok=True)
    index_prefix = index_dir / "ref"

    subprocess.run(
        [
            "bowtie2-build",
            "-f",
            "--threads",
            "1",
            str(ref_fa),
            str(index_prefix),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )

    with TemporaryDirectory(dir=output_dir, prefix="align_") as tmp_raw:
        tmp = Path(tmp_raw)
        sam_file = tmp / "aln.sam"
        subprocess.run(
            [
                "bowtie2",
                "-p",
                "1",
                "--seed",
                str(seed),
                "--no-unal",
                "-x",
                str(index_prefix),
                "-U",
                str(fq_file),
                "-S",
                str(sam_file),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

        cram_file = output_dir / "sample.cram"
        pysam.sort(
            "-O",
            "cram",
            "-o",
            str(cram_file),
            "--reference",
            str(ref_fa),
            "-T",
            str(tmp / "sort"),
            str(sam_file),
        )

    _sanitize_cram_header(cram_file, ref_fa)
    pysam.index(str(cram_file))
    return cram_file


def _make_msp_files(
    output_dir: Path,
    sequences: dict[int, str],
    rng: Generator,
    n_msps: int,
) -> None:
    """Write msp_map.tsv, annotation.tsv, regions.bed and a JSON mirror."""
    msp_names = [f"msp_{i:04d}" for i in range(1, n_msps + 1)]
    gene_ids = list(sequences.keys())
    assignments = rng.choice(msp_names, size=len(gene_ids))

    msp_rows: list[dict[str, str | int]] = []
    annotation_rows: list[dict[str, str | int]] = []
    bed_rows: list[dict[str, int]] = []
    json_map: dict[str, list[int]] = {msp: [] for msp in msp_names}

    for gene_id, msp_name in zip(gene_ids, assignments, strict=True):
        gene_length = len(sequences[gene_id])
        msp_rows.append(
            {
                "msp_name": msp_name,
                "gene_id": gene_id,
                "gene_name": gene_id,
                "gene_category": "core",
            }
        )
        annotation_rows.append(
            {"gene_id": gene_id, "gene_name": gene_id, "gene_length": gene_length}
        )
        bed_rows.append({"gene_id": gene_id, "startpos": 0, "gene_length": gene_length})
        json_map[msp_name].append(gene_id)

    pd.DataFrame(msp_rows).to_csv(
        output_dir / "msp_map.tsv", sep="\t", index=False
    )
    pd.DataFrame(annotation_rows).to_csv(
        output_dir / "annotation.tsv", sep="\t", index=False
    )
    pd.DataFrame(bed_rows).to_csv(
        output_dir / "regions.bed", sep="\t", index=False, header=False
    )

    with (output_dir / "msp_map.json").open("wt", encoding="utf-8") as fh:
        json.dump(json_map, fh, indent=2, sort_keys=True)


def _make_count_file(
    output_dir: Path,
    sequences: dict[int, str],
    gene_counts: dict[int, int],
) -> None:
    """Write the gene-count table variant calling expects."""
    rows = [
        {"gene_id": gene_id, "gene_length": len(seq), "value": gene_counts[gene_id]}
        for gene_id, seq in sequences.items()
    ]
    with lzma.open(output_dir / "sample.tsv.xz", "wt") as fh:
        fh.write("gene_id\tgene_length\tvalue\n")
        for row in rows:
            fh.write(f"{row['gene_id']}\t{row['gene_length']}\t{row['value']}\n")


def _make_reference_json(
    output_dir: Path,
    n_genes: int,
    n_msps: int,
    n_reads: int,
) -> None:
    """Write the *_reference.json meteor consumes."""
    config = {
        "reference_info": {
            "reference_name": "fixture_reference",
            "reference_date": "2024-01-01",
            "database_type": "complete",
        },
        "reference_file": {
            "database_dir": ".",
            "fasta_dir": ".",
            "fasta_filename": "reference.fa.gz",
        },
        "annotation": {
            "gene_id": {
                "filename": "annotation.tsv",
                "version": "1.0",
            },
            "msp": {
                "filename": "msp_map.tsv",
                "version": "1.0",
            },
            "bed": {
                "filename": "regions.bed",
                "version": "1.0",
            },
        },
    }
    with (output_dir / "reference.json").open("wt", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2, sort_keys=True)

    fixture_info = {
        "n_genes": n_genes,
        "n_msps": n_msps,
        "n_reads": n_reads,
    }
    with (output_dir / "fixture_info.json").open("wt", encoding="utf-8") as fh:
        json.dump(fixture_info, fh, indent=2)


def _make_census_stage1(output_dir: Path, n_reads: int) -> None:
    """Write a census_stage_1.json compatible with meteor's counter."""
    config = {
        "sample_info": {
            "sample_name": "sample",
            "sequencing_device": "Illumina",
            "census_status": 0,
            "full_sample_name": "sample",
        },
        "sample_file": {"fastq_file": "sample.fq"},
        "mapping": {
            "mapping_tool": "bowtie2",
            "mapping_tool_version": "2.5.1",
            "mapping_date": "2024-01-01",
            "reference_name": "fixture_reference",
            "trim": 80,
            "alignment_number": 100,
            "mapping_type": "end-to-end",
            "identity_threshold": 0.95,
            "total_read_count": n_reads,
            "mapped_read_count": n_reads,
            "overall_alignment_rate": 100.0,
            "fastq_files": ["sample.fq"],
            "mapping_file": "sample.cram",
        },
    }
    with (output_dir / "sample_census_stage_1.json").open("wt", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2)


@app.command()
def main(
    output_dir: Path = typer.Option(
        Path("tests/data/fixtures"),
        "--output-dir",
        help="Directory where fixture files are written.",
    ),
    seed: int = typer.Option(42, "--seed", help="Random seed for reproducibility."),
    n_genes: int = typer.Option(
        DEFAULT_GENES, "--n-genes", help="Number of reference gene sequences."
    ),
    n_reads: int = typer.Option(
        DEFAULT_READS, "--n-reads", help="Total number of synthetic reads."
    ),
    n_msps: int = typer.Option(
        DEFAULT_MSP_COUNT, "--n-msps", help="Number of MSP groups."
    ),
) -> None:
    """Generate deterministic meteor benchmark fixtures."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = default_rng(seed)

    sequences = _make_reference(output_dir, rng, n_genes)
    gene_counts = _make_reads(output_dir, sequences, rng, n_reads, DEFAULT_READ_LEN)
    _align_reads(output_dir, output_dir / "reference.fa", output_dir / "sample.fq", seed)
    _make_msp_files(output_dir, sequences, rng, n_msps)
    _make_count_file(output_dir, sequences, gene_counts)
    _make_reference_json(output_dir, n_genes, n_msps, n_reads)
    _make_census_stage1(output_dir, n_reads)

    # The FASTQ is an intermediate artifact; keep the working tree tidy.
    (output_dir / "sample.fq").unlink(missing_ok=True)
    shutil.rmtree(output_dir / "bt2_index", ignore_errors=True)

    typer.echo(f"Fixtures written to {output_dir.resolve()}")


if __name__ == "__main__":
    app()
