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
from ..referencebuilder import ReferenceBuilder
from pathlib import Path
from hashlib import md5
import pytest


@pytest.fixture
def builder_defec(tmp_path: Path, datadir: Path) -> ReferenceBuilder:
    meteor = Component
    meteor.ref_dir = tmp_path
    meteor.ref_name = "defect"
    return ReferenceBuilder(meteor, input_fasta=datadir / "defect_catalogue.fasta.gz")


@pytest.fixture
def builder(tmp_path: Path, datadir: Path) -> ReferenceBuilder:
    meteor = Component
    meteor.ref_dir = tmp_path
    meteor.ref_name = "test"
    meteor.threads = 1
    return ReferenceBuilder(meteor, input_fasta=datadir / "test_catalogue.fasta.gz")


def test_read_reference(builder_defec: ReferenceBuilder):
    res = (
        (
            "1",
            80,
            "ATGAAAATGAACCTGCAAAAGGACATGTTTGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT",
        ),
        ("2", 56, "TGCCAGCATGGCGACTTTGTCTTTGATTGGGGCTTTTTTGGGCGGTGGTAGCGCCC"),
        (
            "4",
            74,
            "ATGAAAATGAACCTGCAAAAGGACGATCGCAAATTGCGATACAAGATGTACAAAGATGGTAAAAAGTGGGTGTT",
        ),
    )
    for i, (header, length, sequence) in enumerate(builder_defec.read_reference()):
        assert header == res[i][0]
        assert length == res[i][1]
        assert sequence == res[i][2]


@pytest.mark.parametrize(
    ("annotation_md5", "fasta_md5"),
    (
        pytest.param(
            "be4ea162246d2f23ed8b33bdf9b209d8",
            "55b4a418bd2814f14dd84b7217762b8b",
            id="Accurate output",
        ),
    ),
)
def test_create_reference(
    builder: ReferenceBuilder, annotation_md5: str, fasta_md5: str
):
    builder.create_reference()
    with builder.output_annotation_file.open("rb") as output_annotation:
        assert md5(output_annotation.read()).hexdigest() == annotation_md5
    with builder.output_fasta_file.open("rb") as output_fasta:
        assert md5(output_fasta.read()).hexdigest() == fasta_md5


def test_execute(builder: ReferenceBuilder):
    builder.execute()
    assert builder.output_annotation_file.exists()
    assert builder.output_fasta_file.exists()
    # 6 bt2 file
    assert len(list(builder.output_fasta_file.parent.glob("*.bt2"))) == 6
