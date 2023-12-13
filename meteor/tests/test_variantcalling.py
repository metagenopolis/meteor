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
from configparser import ConfigParser
import pytest


@pytest.fixture
def vc_builder(datadir: Path, tmp_path: Path) -> VariantCalling:
    meteor = Component
    meteor.ref_dir = datadir / "eva71"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.tmp_dir = tmp_path
    ref_ini_file = datadir / "eva71" / "eva71_reference.ini"
    ref_ini = ConfigParser()
    with ref_ini_file.open("rt", encoding="UTF-8") as ref:
        ref_ini.read_file(ref)
    census_ini_file = datadir / "eva71_bench" / "eva71_bench_census_stage_1.ini"
    census_ini = ConfigParser()
    with census_ini_file.open("rt", encoding="UTF-8") as cens:
        census_ini.read_file(cens)
        sample_info = census_ini["sample_info"]
        stage2_dir = tmp_path / sample_info["sample_name"]
        stage2_dir.mkdir(exist_ok=True, parents=True)
        data_dict = {
            "mapped_sample_dir": datadir / "eva71_bench",
            "census": census_ini,
            "directory": stage2_dir,
            "Stage2FileName": stage2_dir / census_ini_file.name,
            "reference": ref_ini,
        }
    return VariantCalling(meteor, data_dict, 100, 0.5)


def test_execute(vc_builder: VariantCalling) -> None:
    vc_builder.execute()
    output_vcf = (
        vc_builder.census["directory"]
        / f"{vc_builder.census['census']['sample_info']['sample_name']}.vcf.gz"
    )
    assert output_vcf.exists()
    output_consensus = (
        vc_builder.census["directory"]
        / f"{vc_builder.census['census']['sample_info']['sample_name']}_consensus.fasta"
    )
    assert output_consensus.exists()
