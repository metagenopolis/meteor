# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Test variant calling"""
from ..session import Component
from ..variantcalling import VariantCalling
from pathlib import Path
import pytest
import json
from pysam import AlignmentFile, FastaFile
import pandas as pd
from pandas.testing import assert_frame_equal
from hashlib import md5


@pytest.fixture(name="vc_builder")
def fixture_vc_builder(datadir: Path, tmp_path: Path) -> VariantCalling:
    meteor = Component
    meteor.ref_dir = datadir / "eva71"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_path
    ref_json_file = datadir / "eva71" / "eva71_reference.json"
    ref_json = {}
    with ref_json_file.open("rt", encoding="UTF-8") as ref:
        ref_json = json.load(ref)
    census_json_file = datadir / "eva71_bench" / "eva71_bench_census_stage_1.json"
    census_json = {}
    with census_json_file.open("rt", encoding="UTF-8") as cens:
        census_json = json.load(cens)
        sample_info = census_json["sample_info"]
        stage3_dir = tmp_path / sample_info["sample_name"]
        stage3_dir.mkdir(exist_ok=True, parents=True)
        data_dict = {
            "mapped_sample_dir": datadir / "eva71_bench",
            "census": census_json,
            "directory": stage3_dir,
            "Stage3FileName": stage3_dir / census_json_file.name,
            "reference": ref_json,
        }
    return VariantCalling(meteor, data_dict, 100, 3, 1, 0.01, 1, 100)


def test_group_consecutive_positions(vc_builder: VariantCalling, datadir) -> None:
    results_df = pd.read_csv(datadir / "expected_output" / "coverage_pos.tsv", sep="\t")
    reads_dict = dict(zip(results_df["Position"], results_df["Count"]))
    expected_output = pd.read_table(
        datadir / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
    )
    df = vc_builder.group_consecutive_positions(reads_dict, "1", 7408)
    df = df.astype(
        {"gene_id": "str", "startpos": "int64", "endpos": "int64", "coverage": "int64"}
    )
    expected_output = expected_output.astype(
        {"gene_id": "str", "startpos": "int64", "endpos": "int64", "coverage": "int64"}
    )
    assert df.equals(expected_output)


def test_count_reads_in_gene(vc_builder: VariantCalling, datadir) -> None:
    cram_file = datadir / "eva71_bench" / "eva71_bench.cram"
    reference_file = datadir / "eva71" / "fasta" / "eva71.fasta.gz"
    with AlignmentFile(str(cram_file.resolve()), "rc") as cram:
        with FastaFile(filename=str(reference_file.resolve())) as Fasta:
            reads_dict = vc_builder.count_reads_in_gene(cram, "1", 7408, Fasta)
            result = pd.DataFrame(reads_dict.items(), columns=["Position", "Count"])
    expected_output = pd.read_table(datadir / "expected_output" / "coverage_pos.tsv")
    assert result.equals(expected_output)


def test_filter_low_cov_sites(vc_builder: VariantCalling, datadir) -> None:
    vc_builder.matrix_file = datadir / "eva71_bench" / "eva71_bench.tsv.xz"
    cram_file = datadir / "eva71_bench" / "eva71_bench.cram"
    reference_file = datadir / "eva71" / "fasta" / "eva71.fasta.gz"

    result_df, _ = vc_builder.filter_low_cov_sites(cram_file, reference_file)
    expected_output = pd.read_table(
        datadir / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
    ).set_index("gene_id")
    assert_frame_equal(
        result_df.sort_values(by=list(result_df.columns)).reset_index(drop=True),
        expected_output.sort_values(by=list(expected_output.columns)).reset_index(
            drop=True
        ),
        check_like=True,
    )


def test_create_consensus(
    vc_builder: VariantCalling, datadir: Path, tmp_path: Path
) -> None:
    vc_builder.matrix_file = datadir / "eva71_bench" / "eva71_bench.tsv.xz"
    vc_builder.meteor.DEFAULT_GAP_CHAR = "?"
    reference_file = datadir / "eva71" / "fasta" / "eva71.fasta.gz"
    consensus_file = tmp_path / "consensus.fasta.xz"
    bed_file = datadir / "eva71" / "database" / "eva71.bed"
    vcf_file = datadir / "eva71_bench" / "eva71_bench.vcf.gz"
    low_cov_sites = pd.read_table(
        datadir / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
    ).set_index("gene_id")

    vc_builder.create_consensus(
        reference_file,
        consensus_file,
        low_cov_sites,
        pd.DataFrame(),
        vcf_file,
        bed_file,
    )
    assert consensus_file.exists()
    with consensus_file.open("rb") as consensus:
        assert md5(consensus.read()).hexdigest() == "3dcc531550bff705949620224d9950f4"


def test_execute(vc_builder: VariantCalling) -> None:
    vc_builder.execute()
    output_vcf = (
        vc_builder.census["directory"]
        / f"{vc_builder.census['census']['sample_info']['sample_name']}.vcf.gz"
    )
    assert output_vcf.exists()
    output_consensus = (
        vc_builder.census["directory"]
        / f"{vc_builder.census['census']['sample_info']['sample_name']}_consensus.fasta.xz"
    )
    assert output_consensus.exists()


def _synthetic_alignments(tmp_path: Path, n_genes: int = 8, covered: int = 5):
    """FASTA of n_genes random genes (300 bp) and a sorted, indexed CRAM with
    reads on the first `covered` genes (higher depth in the middle), plus the
    matching count table."""
    import lzma
    import random
    from pysam import AlignedSegment, AlignmentHeader, faidx, index, sort

    rng = random.Random(5)
    genes = {str(g): "".join(rng.choice("ACGT") for _ in range(300)) for g in range(1, n_genes + 1)}
    fasta = tmp_path / "ref.fasta"
    fasta.write_text("".join(f">{g}\n{s}\n" for g, s in genes.items()))
    faidx(str(fasta))
    header = AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": g, "LN": 300} for g in genes]}
    )
    unsorted = tmp_path / "unsorted.bam"
    counts = {}
    with AlignmentFile(str(unsorted), "wb", header=header) as out:
        for tid, (gene, seq) in enumerate(genes.items()):
            if tid >= covered:
                continue
            starts = list(range(0, 250, 10)) + list(range(100, 150, 5)) * (tid + 1)
            counts[gene] = len(starts)
            for k, start in enumerate(starts):
                aln = AlignedSegment(header)
                aln.query_name = f"r{gene}_{k}"
                aln.reference_id = tid
                aln.reference_start = start
                aln.cigarstring = "50M"
                aln.query_sequence = seq[start : start + 50]
                aln.query_qualities = [30] * 50
                aln.mapping_quality = 40
                aln.set_tag("NM", 0)
                out.write(aln)
    cram = tmp_path / "sample.cram"
    sort("-O", "cram", "--reference", str(fasta), "-o", str(cram), str(unsorted), catch_stdout=False)
    index(str(cram))
    count_table = tmp_path / "sample.tsv.xz"
    with lzma.open(count_table, "wt") as out:
        out.write("gene_id\tgene_length\tvalue\n")
        for gene in genes:
            out.write(f"{gene}\t300\t{counts.get(gene, 0)}\n")
    return fasta, cram, count_table


def test_runs_below_and_depth_of_absent_gene(tmp_path: Path) -> None:
    import numpy as np
    from ..variantcalling import depth_array, runs_below

    starts, ends, values = runs_below(np.array([0, 0, 5, 5, 1, 9, 9, 0]), 3)
    assert starts.tolist() == [0, 4, 7] and ends.tolist() == [2, 5, 8]
    assert values.tolist() == [0, 1, 0]
    fasta, cram, _ = _synthetic_alignments(tmp_path)
    with AlignmentFile(str(cram), "rc", reference_filename=str(fasta)) as aln:
        with FastaFile(str(fasta)) as ref:
            assert depth_array(aln, "absent", 40, ref, 100).tolist() == [0] * 40
            assert depth_array(aln, "1", 300, ref, 100).max() > 3


def test_low_cov_and_ignore_dicts() -> None:
    import numpy as np
    from ..variantcalling import gene_ignore_to_dict, low_cov_to_dict

    assert low_cov_to_dict({"1": "kept"}) == {"1": "kept"}
    assert low_cov_to_dict(pd.DataFrame(columns=["startpos", "endpos"])) == {}
    table = pd.DataFrame(
        {"startpos": [0, 10, 5], "endpos": [3, 20, 8]}, index=["2", "1", "2"]
    )
    result = low_cov_to_dict(table)
    assert np.array_equal(result["2"], [[0, 3], [5, 8]])
    assert np.array_equal(result["1"], [[10, 20]])
    assert gene_ignore_to_dict({1: 10}) == {1: 10}
    assert gene_ignore_to_dict(pd.DataFrame()) == {}
    assert gene_ignore_to_dict(pd.DataFrame({"gene_length": [7.0]}, index=[3])) == {3: 7}


def test_gene_counts(vc_builder: VariantCalling, datadir: Path) -> None:
    vc_builder.matrix_file = datadir / "eva71_bench" / "eva71_bench.tsv.xz"
    counts = vc_builder.gene_counts()
    assert counts and all(value > 0 for value in counts.values())
    assert vc_builder.gene_counts() is counts


def test_create_balanced_bed_chunks(vc_builder: VariantCalling, tmp_path: Path) -> None:
    merged = pd.DataFrame(
        {
            "msp_name": ["m1", "m1", "m2", "m3", "m4"],
            "gene_id": [1, 2, 3, 4, 5],
            "startpos": 0,
            "gene_length": [100, 200, 300, 400, 500],
        }
    )
    chunks = vc_builder.create_balanced_bed_chunks(merged, {1: 50.0, 3: 40.0, 4: 5.0}, 2, tmp_path)
    contents = [pd.read_csv(c, sep="\t", header=None) for c in chunks]
    genes = sorted(g for c in contents for g in c[0])
    assert genes == [1, 2, 3, 4, 5]
    # heaviest chunk first: m1 (51 + 1) alone, then m2 + m3 + m4
    assert sorted(contents[0][0]) == [1, 2]
    # more chunks than MSPs: one chunk per MSP
    assert len(vc_builder.create_balanced_bed_chunks(merged, {}, 10, tmp_path)) == 4


def test_write_marker_reference_and_alignments(vc_builder: VariantCalling, tmp_path: Path) -> None:
    fasta, cram, _ = _synthetic_alignments(tmp_path)
    Path(f"{cram}.crai").unlink()  # rebuilt by reference_ids_in_cram
    assert vc_builder.reference_ids_in_cram(cram) == {0, 1, 2, 3, 4}
    bam, genes = vc_builder.write_marker_alignments(cram, fasta, {2, 7}, tmp_path)
    assert genes == [1, 2, 3, 4, 5, 7]
    with AlignmentFile(str(bam)) as marker, AlignmentFile(
        str(cram), "rc", reference_filename=str(fasta)
    ) as original:
        assert list(marker.references) == ["1", "2", "3", "4", "5", "7"]
        assert [(r.reference_name, r.reference_start) for r in marker] == [
            (r.reference_name, r.reference_start) for r in original
        ]
    reference = vc_builder.write_marker_reference(fasta, [2, 7], tmp_path)
    with FastaFile(reference) as marker_ref, FastaFile(str(fasta)) as ref:
        assert list(marker_ref.references) == ["2", "7"]
        assert marker_ref.fetch("7") == ref.fetch("7")


def test_filter_low_cov_sites_parallel(vc_builder: VariantCalling, tmp_path: Path) -> None:
    """Gene batches in worker processes give the same runs as one process"""
    fasta, cram, counts = _synthetic_alignments(tmp_path)
    vc_builder.matrix_file = counts
    vc_builder.meteor.threads = 1
    single, ignore_single = vc_builder.filter_low_cov_sites(cram, fasta)
    vc_builder.meteor.threads = 2
    try:
        parallel, ignore_parallel = vc_builder.filter_low_cov_sites(cram, fasta)
    finally:
        vc_builder.meteor.threads = 1
    assert len(single) > 0
    assert_frame_equal(single, parallel)
    assert_frame_equal(ignore_single, ignore_parallel)
    assert sorted(ignore_single.index) == [6, 7, 8]
    # genes restricted to a subset without alignment: no run, no error
    empty, _ = vc_builder.filter_low_cov_sites(cram, fasta, gene_subset={8})
    assert len(empty) == 0


def test_create_consensus_ignored_gene(vc_builder: VariantCalling, datadir: Path, tmp_path: Path) -> None:
    import lzma

    vc_builder.meteor.DEFAULT_GAP_CHAR = "?"
    consensus_file = tmp_path / "consensus.fasta.xz"
    vc_builder.create_consensus(
        datadir / "eva71" / "fasta" / "eva71.fasta.gz",
        consensus_file,
        {},
        {1: 7408},
        datadir / "eva71_bench" / "eva71_bench.vcf.gz",
        datadir / "eva71" / "database" / "eva71.bed",
    )
    with lzma.open(consensus_file, "rt") as consensus:
        assert consensus.read() == ">1\n" + "?" * 7408 + "\n"


def test_variant_calling_with_fake_freebayes(vc_builder: VariantCalling, datadir: Path, fake_freebayes) -> None:
    """execute() end to end, freebayes replaced by a script printing the
    reference VCF of the test data; a second run reuses the VCF and the
    low coverage regions."""
    fake_freebayes.setenv("FAKE_FREEBAYES_VCF", str(datadir / "eva71_bench" / "eva71_bench.vcf.gz"))
    vc_builder.meteor.DEFAULT_GAP_CHAR = "?"
    vc_builder.execute()
    sample = vc_builder.census["census"]["sample_info"]["sample_name"]
    directory = vc_builder.census["directory"]
    consensus = directory / f"{sample}_consensus.fasta.xz"
    for output in (directory / f"{sample}.vcf.gz", directory / f"{sample}.pickle", consensus):
        assert output.exists()
    assert vc_builder.census["Stage3FileName"].exists()
    first = consensus.read_bytes()
    # same file as VariantCalling.execute of meteor 2.0.22 with this VCF
    assert md5(first).hexdigest() == "916fe0018af1eebc59592a75c770570b"
    consensus.unlink()
    vc_builder.execute()
    assert consensus.read_bytes() == first


def test_variant_calling_failed_freebayes(vc_builder: VariantCalling, fake_freebayes) -> None:
    fake_freebayes.setenv("FAKE_FREEBAYES_FAIL", "1")
    with pytest.raises(SystemExit):
        vc_builder.execute()
