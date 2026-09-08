"""Parity tests for the Rust-accelerated consensus builder."""

from __future__ import annotations

import lzma
import tempfile
from pathlib import Path

import meteor_core


def _make_fasta(path: Path, sequences: dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in sequences.items():
            fh.write(f">{name}\n{seq}\n")


def _make_bed(path: Path, rows: list[tuple[str, int, int]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for chrom, start, end in rows:
            fh.write(f"{chrom}\t{start}\t{end}\n")


def _make_vcf(path: Path, header_lines: list[str], records: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for line in header_lines:
            fh.write(line + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for rec in records:
            fh.write(rec + "\n")


def test_create_consensus_basic() -> None:
    seq1 = "ACGTACGTACGTACGTACGT"
    seq2 = "AAAAAAAAAA"
    with tempfile.TemporaryDirectory() as tmp_raw:
        tmp = Path(tmp_raw)
        ref = tmp / "ref.fa"
        bed = tmp / "regions.bed"
        vcf = tmp / "variants.vcf"

        _make_fasta(ref, {"1": seq1, "2": seq2})
        _make_bed(bed, [("1", 0, len(seq1)), ("2", 0, len(seq2))])
        _make_vcf(
            vcf,
            [
                '##contig=<ID=1,length=20>',
                '##contig=<ID=2,length=10>',
                '##INFO=<ID=RO,Number=1,Type=Integer,Description="Ref observations">',
                '##INFO=<ID=AO,Number=A,Type=Integer,Description="Alt observations">',
            ],
            [
                # ref_freq = 1/(1+10) < 0.5 -> keep alt only
                "1\t3\t.\tA\tG\t.\t.\tRO=1;AO=10\tGT\t.",
                # ref_freq = 10/(10+1) >= 0.5 -> keep ref + alt
                "1\t5\t.\tC\tT\t.\t.\tRO=10;AO=1\tGT\t.",
                # ref_freq = 0/(0+5) < 0.5 -> keep alt only
                "2\t2\t.\tA\tC\t.\t.\tRO=0;AO=5\tGT\t.",
            ],
        )

        low_cov_sites = [(1, 10, 12)]  # mask positions 10 and 11 (0-based)
        gene_ignore = []

        results = meteor_core.create_consensus(
            str(ref),
            str(vcf),
            str(bed),
            low_cov_sites,
            gene_ignore,
            min_frequency=0.5,
            gap_char="N",
        )
        results_dict = dict(results)

        expected1 = list(seq1)
        expected1[2] = "G"  # alt at 1:3 -> 0-based 2
        expected1[4] = "Y"  # ref+alt at 1:5 -> 0-based 4 -> IUPAC Y
        expected1[10] = "N"
        expected1[11] = "N"

        expected2 = list(seq2)
        expected2[1] = "C"  # alt at 2:2 -> 0-based 1

        assert results_dict[1] == "".join(expected1)
        assert results_dict[2] == "".join(expected2)


def test_create_consensus_ignore_and_low_cov() -> None:
    seq1 = "ACGTACGTACGT"
    with tempfile.TemporaryDirectory() as tmp_raw:
        tmp = Path(tmp_raw)
        ref = tmp / "ref.fa"
        bed = tmp / "regions.bed"
        vcf = tmp / "variants.vcf"

        _make_fasta(ref, {"1": seq1})
        _make_bed(bed, [("1", 0, len(seq1))])
        _make_vcf(
            vcf,
            [
                '##contig=<ID=1,length=12>',
                '##INFO=<ID=RO,Number=1,Type=Integer,Description="Ref observations">',
                '##INFO=<ID=AO,Number=A,Type=Integer,Description="Alt observations">',
            ],
            [],
        )

        results = meteor_core.create_consensus(
            str(ref),
            str(vcf),
            str(bed),
            low_cov_sites=[],
            gene_ignore=[(1, len(seq1))],
            min_frequency=0.5,
            gap_char="N",
        )
        assert dict(results)[1] == "N" * len(seq1)
