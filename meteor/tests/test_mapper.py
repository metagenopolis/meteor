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

"""Test mapper builder main object"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..mapper import Mapper
from pathlib import Path
import pytest
import json


@pytest.fixture
def mapping_builder(datadir: Path, tmp_path: Path) -> Mapper:
    meteor = Component
    meteor.ref_dir = datadir / "eva71"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_path
    ref_json_file = datadir / "eva71" / "eva71_reference.json"
    ref_json = {}
    with ref_json_file.open("rt", encoding="UTF-8") as ref:
        ref_json = json.load(ref)
    census_ini_file = datadir / "eva71_bench_census_stage_0.json"
    census_ini = {}
    with census_ini_file.open("rt", encoding="UTF-8") as cens:
        census_ini = json.load(cens)
        sample_info = census_ini["sample_info"]
        stage1_dir = tmp_path / sample_info["sample_name"]
        stage1_dir.mkdir(exist_ok=True, parents=True)
        data_dict = {
            "census": census_ini,
            "directory": stage1_dir,
            "Stage1FileName": stage1_dir
            / census_ini_file.name.replace("stage_0", "stage_1"),
            "reference": ref_json,
        }
    return Mapper(
        meteor,
        data_dict,
        [str(datadir / "eva71_bench.fq.gz")],
        "end-to-end",
        80,
        10000,
    )


# def test_create_cram(mapping_builder: Mapper,  datadir: Path, tmp_path: Path) -> None:
#     input_sam = datadir / "eva71_bench.sam"
#     output_cram = tmp_path / "eva71_bench.cram"
#     mapping_builder.create_cram(str(input_sam.resolve()),  str(output_cram.resolve()))
#     assert output_cram.exists()


# def test_sort_cram(mapping_builder: Mapper,  datadir: Path, tmp_path: Path) -> None:
#     input_cram = datadir / "eva71_bench.cram"
#     output_cram = tmp_path / "eva71_bench_sorted.cram"
#     mapping_builder.sort_cram(str(input_cram.resolve()), str(output_cram.resolve()))
#     assert output_cram.exists()


def test_execute(mapping_builder: Mapper) -> None:
    mapping_builder.execute()
    output_cram = (
        mapping_builder.census["directory"]
        / f"{mapping_builder.census['census']['sample_info']['sample_name']}_raw.cram"
    )
    assert output_cram.exists()
