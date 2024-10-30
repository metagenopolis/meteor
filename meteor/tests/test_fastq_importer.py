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

"""Import and prepare fastq"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..fastqimporter import FastqImporter
from pathlib import Path
import pytest


@pytest.fixture
def builder(tmp_path: Path, datadir: Path) -> FastqImporter:
    meteor = Component
    meteor.fastq_dir = tmp_path
    return FastqImporter(meteor, datadir, False, "\\d+")


@pytest.fixture
def builder_single(tmp_path: Path, datadir: Path) -> FastqImporter:
    meteor = Component
    meteor.fastq_dir = tmp_path / "single"
    return FastqImporter(meteor, datadir / "single", False, None)


@pytest.fixture
def builder_paired(tmp_path: Path, datadir: Path) -> FastqImporter:
    meteor = Component
    meteor.fastq_dir = tmp_path / "paired"
    return FastqImporter(meteor, datadir / "paired", True, None)


@pytest.fixture
def builder_single_mask(tmp_path: Path, datadir: Path) -> FastqImporter:
    meteor = Component
    meteor.fastq_dir = tmp_path / "single_mask"
    return FastqImporter(meteor, datadir / "single_mask", False, "Zymo_\\d+")


@pytest.fixture
def builder_paired_mask(tmp_path: Path, datadir: Path) -> FastqImporter:
    meteor = Component
    meteor.fastq_dir = tmp_path / "paired_mask"
    return FastqImporter(meteor, datadir / "paired_mask", True, "Zymo_\\d+")


@pytest.fixture
def expected_dir() -> tuple[Path, Path, Path]:
    return Path("test"), Path("test2"), Path("test3")


@pytest.fixture
def expected_dir_mask() -> tuple[Path, Path]:
    return Path("Zymo_6300"), Path("Zymo_6400")


@pytest.mark.parametrize(
    ("fastq_filename", "name"),
    (
        ("test.fastq.gz", "test"),
        ("simple_case.fastq", "simple_case"),
        ("pretty.complex_pain.fq.xz", "pretty.complex_pain"),
        pytest.param(
            "pretty.complex_pain.fasta", "pretty.complex_pain.fasta", id="fasta"
        ),
    ),
)
def test_replace_ext(builder: FastqImporter, fastq_filename: str, name: str) -> None:
    assert builder.replace_ext(fastq_filename) == name


@pytest.mark.parametrize(
    ("fastq_filename"),
    (
        ("test.fastq.gz"),
        pytest.param("pretty.complex_pain.fasta", id="fasta"),
    ),
)
def test_get_tag_none(builder: FastqImporter, fastq_filename: str) -> None:
    assert builder.get_tag(fastq_filename) is None


@pytest.mark.parametrize(
    ("fastq_filename", "tag"),
    (("simple_case_R1.fastq", "1"), ("pretty.complex_pain_2.fq.xz", "2")),
)
def test_get_tag(builder: FastqImporter, fastq_filename: str, tag: str) -> None:
    assert builder.get_tag(fastq_filename) == tag


@pytest.mark.parametrize(
    ("fastq_filename", "tag", "sample_name"),
    (
        ("simple_case_R1.fastq", "1", "simple_case"),
        ("pretty.complex_pain_2.fq.xz", "2", "pretty.complex_pain"),
    ),
)
def test_get_paired_dirname(
    builder: FastqImporter, fastq_filename: str, tag: str, sample_name: str
) -> None:
    assert builder.get_paired_dirname(fastq_filename, tag) == sample_name


def test_execute_single(builder_single: FastqImporter, expected_dir: tuple) -> None:
    builder_single.execute()
    assert all(
        Path(builder_single.meteor.fastq_dir / dir).exists() for dir in expected_dir
    )
    assert all(
        Path(builder_single.meteor.fastq_dir / dir).is_dir() for dir in expected_dir
    )
    assert all(
        len(list(Path(builder_single.meteor.fastq_dir / dir).glob("*.f*q*"))) == 1
        for dir in expected_dir
    )


def test_execute_paired(builder_paired: FastqImporter, expected_dir: tuple) -> None:
    builder_paired.execute()
    assert all(
        Path(builder_paired.meteor.fastq_dir / dir).exists() for dir in expected_dir
    )
    assert all(
        Path(builder_paired.meteor.fastq_dir / dir).is_dir() for dir in expected_dir
    )
    assert all(
        len(list(Path(builder_paired.meteor.fastq_dir / dir).glob("*.f*q*"))) == 2
        for dir in expected_dir
    )


def test_execute_single_mask(
    builder_single_mask: FastqImporter, expected_dir_mask: tuple
) -> None:
    builder_single_mask.execute()
    assert all(
        Path(builder_single_mask.meteor.fastq_dir / dir).exists()
        for dir in expected_dir_mask
    )
    assert all(
        Path(builder_single_mask.meteor.fastq_dir / dir).is_dir()
        for dir in expected_dir_mask
    )
    assert all(
        len(list(Path(builder_single_mask.meteor.fastq_dir / dir).glob("*.f*q*"))) == 1
        for dir in expected_dir_mask
    )


def test_execute_paired_mask(
    builder_paired_mask: FastqImporter, expected_dir_mask: tuple
) -> None:
    builder_paired_mask.execute()
    for dir in expected_dir_mask:
        print(builder_paired_mask.meteor.fastq_dir / dir)
    assert all(
        Path(builder_paired_mask.meteor.fastq_dir / dir).exists()
        for dir in expected_dir_mask
    )
    assert all(
        Path(builder_paired_mask.meteor.fastq_dir / dir).is_dir()
        for dir in expected_dir_mask
    )
    assert all(
        len(list(Path(builder_paired_mask.meteor.fastq_dir / dir).glob("*.f*q*"))) == 2
        for dir in expected_dir_mask
    )
