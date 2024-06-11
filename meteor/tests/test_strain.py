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

"""Test strain analysis"""

from ..session import Component
from ..strain import Strain
from pathlib import Path
import pytest
import pandas as pd


@pytest.fixture
def strain_builder(datadir: Path, tmp_path: Path) -> Strain:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.strain_dir = datadir / "strain"
    meteor.ref_dir = datadir / "catalogue" / "mock"
    meteor.ref_name = "test"
    meteor.threads = 1
    # meteor.mapping_dir = tmp_path
    return Strain(meteor, 100, 3, 3, 0.5, 3, 0.8, True)


def test_filter_coverage(
    strain_builder: Strain,
    datadir: Path,
) -> None:
    filtered_cov = strain_builder.filter_coverage(
        strain_builder.meteor.strain_dir / "test.cram",
        strain_builder.meteor.ref_dir / "database" / "mock.bed",
        strain_builder.meteor.ref_dir / "fasta" / "mock.fasta.gz",
    )
    expected_output = pd.read_table(datadir / "expected_output" / "filtered_cov.tsv")
    assert filtered_cov.reset_index(drop=True).equals(expected_output)
