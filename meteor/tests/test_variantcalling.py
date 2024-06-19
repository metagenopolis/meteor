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
from pysam import AlignmentFile
import pandas as pd
from subprocess import run


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
    return VariantCalling(meteor, data_dict, 100, 3, 3, 0.5)


def test_count_reads_in_gene(vc_builder: VariantCalling, datadir) -> None:
    cram_file = datadir / "eva71_bench" / "eva71_bench.cram"
    reference_file = datadir / "eva71" / "fasta" / "eva71.fasta.gz"
    with AlignmentFile(str(cram_file.resolve()), "rc") as cram:
        reads_dict = vc_builder.count_reads_in_gene(cram, "1", 7408, reference_file)
        result = pd.DataFrame(reads_dict.items(), columns=["Position", "Count"])
    expected_output = pd.read_table(datadir / "expected_output" / "coverage_pos.tsv")
    assert result.equals(expected_output)


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


def test_filter_low_cov_sites_python(
    vc_builder: VariantCalling, datadir, tmpdir
) -> None:
    vc_builder.matrix_file = datadir / "eva71_bench" / "eva71_bench.tsv.xz"
    cram_file = datadir / "eva71_bench" / "eva71_bench.cram"
    reference_file = datadir / "eva71" / "fasta" / "eva71.fasta.gz"
    result = tmpdir / "filt_cov.tsv"
    with result.open("wt") as out:
        vc_builder.filter_low_cov_sites_python(cram_file, reference_file, out)
    result_df = pd.read_table(
        result, names=["gene_id", "startpos", "endpos", "coverage"], sep="\t"
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
    )
    assert result_df.equals(expected_output)


def test_filter_low_cov_sites(vc_builder: VariantCalling, datadir, tmpdir) -> None:
    vc_builder.matrix_file = datadir / "eva71_bench" / "eva71_bench.tsv.xz"
    cram_file = datadir / "eva71_bench" / "eva71_bench.cram"
    result = tmpdir / "filt_cov.tsv"
    output = run(
        [
            "bedtools",
            "genomecov",
            "-bga",
            "-ibam",
            str(cram_file.resolve()),
        ],
        check=False,
        capture_output=True,
    ).stdout.decode("utf-8")
    with result.open("wt") as out:
        vc_builder.filter_low_cov_sites(output, out)
    result_df = pd.read_table(
        result, names=["gene_id", "startpos", "endpos", "coverage"], sep="\t"
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "coverage_expected.tsv", header=0, sep="\t"
    )
    assert result_df.equals(expected_output)


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
