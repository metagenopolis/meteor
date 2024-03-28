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
import pandas as pd


@pytest.fixture
def vc_builder(datadir: Path, tmp_path: Path) -> VariantCalling:
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
    return VariantCalling(meteor, data_dict, 100, 3, 0.5)


# def test_get_regions(vc_builder: VariantCalling) -> None:
#     regions = vc_builder.get_regions()
#     print(regions)


# def test_count_reads_in_gene(vc_builder: VariantCalling) -> None:
#     reads_in_gene = vc_builder.count_reads_in_gene("ENSG00000279457")


# def test_filter_low_cov_sites(vc_builder: VariantCalling, datadir: Path) -> None:
#     vc_builder.filter_low_cov_sites(cram_file)
#     expected_output = pd.read_table(datadir / "expected_output" / "D15B1_low_cov.tsv")
#     assert profiler_standard.gene_count.equals(expected_output)


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
