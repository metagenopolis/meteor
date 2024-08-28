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
from pathlib import Path
from hashlib import md5
from urllib.request import urlretrieve
from time import time
import tarfile
import json
from typing import ClassVar
import sys


@dataclass
class Downloader(Session):
    """Download and prepare catalogues"""

    CONFIG_DATA_FILE: ClassVar[Path] = Path("data/zenodo.json")
    TEST_CATALOGUE: ClassVar[str] = "test"

    meteor: type[Component]
    choice: str
    taxonomy: bool
    check_md5: bool
    data_type: str = "file_info"
    catalogues_config: dict = field(default_factory=dict)
    start_time: float = field(default_factory=float)

    @staticmethod
    def load_catalogues_config() -> dict:
        try:
            config_data = importlib.resources.files("meteor") / str(
                Downloader.CONFIG_DATA_FILE
            )
            with importlib.resources.as_file(config_data) as configuration_path:
                with configuration_path.open("rt", encoding="UTF-8") as config:
                    return json.load(config)
        except FileNotFoundError:
            logging.error(
                "The file %s is missing in meteor source",
                Downloader.CONFIG_DATA_FILE.name,
            )
            sys.exit(1)

    @staticmethod
    def get_available_catalogues() -> list[str]:
        catalogues_config = Downloader.load_catalogues_config()
        available_catalogues = list(catalogues_config.keys())
        available_catalogues.remove(Downloader.TEST_CATALOGUE)
        return available_catalogues

    def __post_init__(self) -> None:
        self.catalogues_config = Downloader.load_catalogues_config()
        self.meteor.ref_dir.mkdir(exist_ok=True, parents=True)
        if self.taxonomy:
            self.data_type = "taxonomy_info"

    def getmd5(self, catalog: Path) -> str:
        """Compute in md5 chunck by chunk to avoid memory overload

        :param catalog: (Path) A path object to the downloaded catalog
        :return: (str) A string of the md5sum
        """
        logging.info("Checking md5 %s", str(catalog))
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

    def execute(self) -> None:
        try:
            # for choice in self.user_choice:
            logging.info(
                "Download %s catalogue",
                self.catalogues_config[self.choice][self.data_type]["filename"],
            )
            url = self.catalogues_config[self.choice][self.data_type]["catalogue"]
            md5_expect = self.catalogues_config[self.choice][self.data_type]["md5"]
            catalogue = (
                self.meteor.ref_dir
                / self.catalogues_config[self.choice][self.data_type]["filename"]
            )
            urlretrieve(url, filename=catalogue, reporthook=self.show_progress)
            print(flush=True)
            if self.choice == Downloader.TEST_CATALOGUE:
                for sample in self.catalogues_config[self.choice]["samples"]:
                    logging.info("Download %s fastq file", sample)
                    url_fastq = self.catalogues_config[self.choice]["samples"][sample][
                        "catalogue"
                    ]
                    fastq_test = (
                        self.meteor.tmp_dir
                        / self.catalogues_config[self.choice]["samples"][sample][
                            "filename"
                        ]
                    )
                    md5fastq_expect = self.catalogues_config[self.choice]["samples"][
                        sample
                    ]["md5"]
                    urlretrieve(
                        url_fastq, filename=fastq_test, reporthook=self.show_progress
                    )
                    print(flush=True)
                    assert md5fastq_expect == self.getmd5(fastq_test)
            if self.check_md5:
                assert md5_expect == self.getmd5(catalogue)
            self.extract_tar(catalogue)
            logging.info(
                "The catalogue is now ready to be used in the folder: %s",
                str(catalogue.with_suffix("").stem),
            )
        except AssertionError:
            logging.error("MD5sum of %s has a different value than expected", catalogue)
            sys.exit(1)
