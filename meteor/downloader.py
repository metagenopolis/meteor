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

"""Download and index reference"""

import logging
import importlib.resources
from dataclasses import dataclass, field
from meteor.session import Session, Component
from configparser import ConfigParser
from pathlib import Path
from hashlib import md5
from urllib.request import urlretrieve
from typing import Type
from time import time
import tarfile


@dataclass
class Downloader(Session):
    """Download and prepare catalogues"""

    meteor: Type[Component]
    choice: str
    check_md5: bool
    catalogues_config: ConfigParser = field(default_factory=ConfigParser)
    start_time: float = field(default_factory=float)

    def __post_init__(self) -> None:
        try:
            config_data = (
                importlib.resources.files("meteor") / "data/dataverse_inrae.ini"
            )
            with importlib.resources.as_file(config_data) as config:
                self.catalogues_config.read_file(config)
        except AssertionError:
            logging.error("The file dataverse_inrae.ini is missing in meteor source")
        self.meteor.ref_dir.mkdir(exist_ok=True)

    def getmd5(self, catalog: Path) -> str:
        """Compute in md5 chunck by chunk to avoid memory overload

        :param catalog: (Path) A path object to the downloaded catalog
        :return: (str) A string of the md5sum
        """
        logging.info("Checking file md5 %s", str(catalog))
        with catalog.open("rb") as f:
            file_hash = md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)
        return file_hash.hexdigest()

    def show_progress(
        self, block_num: int, block_size: int, total_size: int
    ) -> None:  # pragma: no cover
        """Show download progress block per block

        :param block_num: Number of the block
        :param block_size: One block size
        :param total_size: Total size of the file in block
        """
        if block_num == 0:
            self.start_time = time()
            return
        duration = time() - self.start_time
        progress_size = block_num * block_size
        speed = progress_size / (1024 * duration)
        percent = block_num * block_size / total_size
        hours, rem = divmod(duration, 3600)
        minutes, seconds = divmod(rem, 60)
        print(
            f"Download of {self.choice} catalogue : {percent:.1%}, {progress_size/ (1024 * 1024): 8.0f} MB, "
            f"{speed: 6.0f} KB/s, {int(hours):0>2}:{int(minutes):0>2}:{int(seconds):0>2} elapsed time.",
            end="\r",
            flush=True,
        )

    def extract_tar(self, catalogue: Path) -> None:
        """Extract tar file

        :param catalogue: (Path) A path object to the given catalog
        """
        logging.info("Extracting %s catalogue", self.choice)
        with tarfile.open(catalogue) as tar:
            tar.extractall(path=self.meteor.ref_dir)
        catalogue.unlink(missing_ok=True)

    def execute(self) -> bool:
        try:
            # for choice in self.user_choice:
            logging.info("Download %s microbiome reference catalogue", self.choice)
            url = self.catalogues_config[self.choice]["catalogue"]
            md5_expect = self.catalogues_config[self.choice]["md5"]
            catalogue = (
                self.meteor.ref_dir / self.catalogues_config[self.choice]["filename"]
            )
            urlretrieve(url, filename=catalogue, reporthook=self.show_progress)
            print(flush=True)
            if self.choice == "test":
                logging.info("Download test fastq file")
                url_fastq = self.catalogues_config[self.choice]["fastq"]
                fastq_test = (
                    self.meteor.tmp_dir
                    / self.catalogues_config[self.choice]["fastqfilename"]
                )
                md5fastq_expect = self.catalogues_config[self.choice]["md5fastq"]
                urlretrieve(
                    url_fastq, filename=fastq_test, reporthook=self.show_progress
                )
                print(flush=True)
                assert md5fastq_expect == self.getmd5(fastq_test)
            if self.check_md5:
                assert md5_expect == self.getmd5(catalogue)
            self.extract_tar(catalogue)
            logging.info("Catalogue %s is now ready to be used.", self.choice)
        except AssertionError:
            logging.error("MD5sum of %s has a different value than expected", catalogue)
        return True
