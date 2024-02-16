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
from ..downloader import Downloader
from pathlib import Path
import pytest


@pytest.fixture
def downer(tmp_path: Path) -> Downloader:
    meteor = Component
    meteor.fastq_dir = tmp_path
    meteor.ref_dir = tmp_path
    meteor.tmp_dir = tmp_path
    meteor.ref_name = "test"
    meteor.threads = 1
    return Downloader(meteor, choice="test", taxonomy=False, check_md5=True)


def test_check_md5(downer: Downloader, datadir: Path) -> None:
    catalogue_file = datadir / "test.tar.xz"
    assert downer.getmd5(catalogue_file) == "ccd7327b5badc2b612452851e8a67ec1"


def test_extract_tar(downer: Downloader, datadir: Path) -> None:
    catalogue_file = datadir / "test.tar.xz"
    downer.extract_tar(catalogue_file)
    assert Path(downer.meteor.ref_dir / "test").exists()


def test_execute(downer: Downloader) -> None:
    downer.execute()
    assert Path(downer.meteor.ref_dir / "mock").exists()
    assert Path(downer.meteor.ref_dir / "test.fastq.gz").exists()
