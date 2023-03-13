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

"""Test reference builder main objects"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..mapper2 import Mapper2
from pathlib import Path
from hashlib import md5
import pytest


@pytest.fixture
def mapping_builder(tmp_path: Path, datadir: Path):
    meteor = Component
    meteor.ref_dir = tmp_path
    meteor.ref_name = "test"
    meteor.threads = 1
    return Mapper2(meteor, {}, [str(datadir / "eva71.fq.xz")], "end-to-end", 80, 10000, "smart_shared_reads")


def test_create_bam(mapping_builder: Mapper2,  datadir: Path, tmp_path: Path):
    input_sam = datadir / "eva71.sam"
    result_bam = tmp_path / "eva71.bam"
    mapping_builder.create_bam(str(input_sam.resolve()),  str(result_bam.resolve()))
    with result_bam.open("rb") as bam:
        assert  md5(bam.read()).hexdigest() == "427347d53952911c850bd95e577175a8"
