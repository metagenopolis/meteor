"""Parity tests for Rust-accelerated CRAM coverage pileup."""

from __future__ import annotations

from pathlib import Path

import pytest
from pysam import AlignmentFile, FastaFile

meteor_core = pytest.importorskip("meteor_core")

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "tests" / "data" / "fixtures"
CRAM = FIXTURE_DIR / "sample.cram"
REF = FIXTURE_DIR / "reference.fa"
MAX_DEPTH = 8000

GENES = ["1", "2", "3", "10", "42", "100"]


def _pysam_counts(gene_name: str, gene_length: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    with AlignmentFile(
        str(CRAM), "rc", reference_filename=str(REF), threads=1
    ) as cram, FastaFile(str(REF)) as fasta:
        for column in cram.pileup(
            contig=gene_name,
            start=0,
            end=gene_length,
            stepper="all",
            max_depth=MAX_DEPTH,
            fastafile=fasta,
            multiple_iterators=False,
        ):
            counts[column.reference_pos] = sum(
                1 for read in column.pileups if not read.is_del and not read.is_refskip
            )
    return counts


def _gene_length(gene_name: str) -> int:
    with FastaFile(str(REF)) as fasta:
        return fasta.get_reference_length(gene_name)


@pytest.mark.parametrize("gene_name", GENES)
def test_count_reads_in_gene_parity(gene_name: str) -> None:
    gene_length = _gene_length(gene_name)
    expected = _pysam_counts(gene_name, gene_length)
    rust_counts = meteor_core.count_reads_in_gene(
        str(CRAM), str(REF), gene_name, gene_length, MAX_DEPTH
    )
    assert rust_counts == expected
