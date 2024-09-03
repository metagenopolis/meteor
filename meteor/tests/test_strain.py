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
from hashlib import md5
import pytest
import pandas as pd


@pytest.fixture
def strain_builder(datadir: Path, tmp_path: Path) -> Strain:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.strain_dir = tmp_path / "strain"
    meteor.ref_dir = datadir / "catalogue" / "mock"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapped_sample_dir = datadir / "map" / "test"
    meteor.DEFAULT_GAP_CHAR = "?"
    json_data = {
        "directory": tmp_path,
        "mapped_sample_dir": datadir / "map" / "test",
        "census": {"sample_info": {"sample_name": "test"}},
    }
    # meteor.mapping_dir = tmp_path
    # return Strain(meteor, 100, 3, 3, 0.5, 1, 3, 0.8, True, False, json_data=json_data)
    return Strain(meteor, 100, 3, 3, 0.5, 1, 1, 0.8, 100, True, json_data=json_data)


def test_filter_coverage(
    strain_builder: Strain,
    datadir: Path,
) -> None:
    filtered_cov = strain_builder.filter_coverage(
        strain_builder.json_data["mapped_sample_dir"] / "test.cram",
        strain_builder.meteor.ref_dir / "database" / "mock.bed",
        strain_builder.meteor.ref_dir / "fasta" / "mock.fasta.gz",
    )
    expected_output = pd.read_table(datadir / "expected_output" / "filtered_cov.tsv")
    assert filtered_cov.reset_index(drop=True).equals(expected_output)


@pytest.mark.parametrize(
    "sequence,expected",
    [
        ("?????", True),
        ("?", True),
        ("?? ??", False),  # contains space
        (
            "",
            False,
        ),
        ("abc?", False),
        ("??a?", False),
        ("a????", False),
    ],
)
def test_is_only_question_marks(strain_builder: Strain, sequence, expected):
    assert strain_builder.is_only_question_marks(sequence) == expected


def test_get_msp_variant(strain_builder, datadir: Path, tmp_path: Path):
    consensus_file = datadir / "strain" / "test_consensus.fasta.xz"
    msp_file = strain_builder.meteor.ref_dir / "database" / "mock_genomes.tsv"
    cram_file = strain_builder.json_data["mapped_sample_dir"] / "test.cram"
    bed_file = strain_builder.meteor.ref_dir / "database" / "mock.bed"
    reference_file = strain_builder.meteor.ref_dir / "fasta" / "mock.fasta.gz"
    strain_builder.get_msp_variant(
        consensus_file, msp_file, cram_file, bed_file, reference_file
    )
    BS = tmp_path / "BS.fasta.xz"
    PA = tmp_path / "PA.fasta.xz"
    assert BS.exists()
    assert PA.exists()
    with BS.open("rb") as out:
        # assert md5(out.read()).hexdigest() == "958199bf4afa9adb56c9e0100a0cc23c"
        assert md5(out.read()).hexdigest() == "c0f943768ca8629b8be4a995bb3af341"
    with PA.open("rb") as out:
        # assert md5(out.read()).hexdigest() == "cec667aad0449853b3bd7db689dd7ed3"
        assert md5(out.read()).hexdigest() == "cec667aad0449853b3bd7db689dd7ed3"


def test_execute(strain_builder, tmp_path: Path) -> None:
    strain_builder.execute()
    BS = tmp_path / "strain" / "test" / "BS.fasta.xz"
    assert BS.exists()
    with BS.open("rb") as out:
        assert md5(out.read()).hexdigest() == "41d0686c1dd663d5511666e4b169d33c"
