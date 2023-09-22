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
from ..mapper2 import Mapper2
from pathlib import Path
from configparser import ConfigParser
import pytest


@pytest.fixture
def mapping_builder(datadir: Path, tmp_path: Path) -> Mapper2:
    meteor = Component
    meteor.ref_dir = datadir / "eva71"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_path
    ref_ini_file = datadir / "eva71" / "eva71_reference.ini"
    ref_ini = ConfigParser()
    with ref_ini_file.open("rt", encoding="UTF-8") as ref:
        ref_ini.read_file(ref)
    census_ini_file = datadir / "eva71_bench_census_stage_0.ini"
    census_ini = ConfigParser()
    with census_ini_file.open("rt", encoding="UTF-8") as cens:
        census_ini.read_file(cens)
        sample_info = census_ini["sample_info"]
        stage1_dir = tmp_path / sample_info["sample_name"]
        stage1_dir.mkdir(exist_ok=True, parents=True)
        data_dict = {
            "census": census_ini,
            "directory": stage1_dir,
            "Stage1FileName": stage1_dir / census_ini_file.name.replace("stage_0", "stage_1"),
            "reference": ref_ini
        }
    return Mapper2(meteor, data_dict, [str(datadir / "eva71_bench.fq.gz")],
                   "end-to-end", 80, 10000, "smart_shared_reads")


# def test_create_bam(mapping_builder: Mapper2,  datadir: Path, tmp_path: Path) -> None:
#     input_sam = datadir / "eva71_bench.sam"
#     output_bam = tmp_path / "eva71_bench.bam"
#     mapping_builder.create_bam(str(input_sam.resolve()),  str(output_bam.resolve()))
#     assert output_bam.exists()


# def test_sort_bam(mapping_builder: Mapper2,  datadir: Path, tmp_path: Path) -> None:
#     input_bam = datadir / "eva71_bench.bam"
#     output_bam = tmp_path / "eva71_bench_sorted.bam"
#     mapping_builder.sort_bam(str(input_bam.resolve()), str(output_bam.resolve()))
#     assert output_bam.exists()


def test_execute(mapping_builder: Mapper2) -> None:
    mapping_builder.execute()
    output_sam = (mapping_builder.census["directory"] /
                  f"{mapping_builder.census['census']['sample_info']['sample_name']}.sam")
    assert output_sam.exists()
