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

from dataclasses import dataclass, field
from session import Session, Component
from referencebuilder import ReferenceBuilder
from configparser import ConfigParser
from pathlib import Path
from hashlib import md5
from urllib.request import urlretrieve
from typing import Type
import logging

@dataclass
class Downloader(Session):
    """Download and prepare catalogues"""
    meteor: Type[Component]
    choice: str
    configuration_path: Path = field(default_factory=Path)
    catalogues_config: ConfigParser = field(default_factory=ConfigParser)

    def __post_init__(self)->None:
        try:
            self.configuration_path = Path(__file__).parent / "dataverse_inrae.ini"
            assert self.configuration_path.exists()
        except AssertionError:
            logging.error("The file dataverse_inrae.ini is missing in meteor source")
        with self.configuration_path.open("rt", encoding="UTF-8") as config:
            self.catalogues_config.read_file(config)
        self.meteor.ref_dir.mkdir(exist_ok=True)

    def getmd5(self, catalog: Path)->str:
        """Compute in md5 chunck by chunk to avoid memory overload
        :param catalog: A path object to the downloaded catalog
        :return: A string of the md5sum
        """
        logging.info("Compute md5sum")
        with catalog.open("rb") as f:
            file_hash = md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)
        return file_hash.hexdigest()

    def show_progress(self, block_num:int, block_size:int, total_size:int):
        """Show download progress block per block"""
        print(f"Download of {self.choice} catalogue : {round(block_num * block_size / total_size *100,2)}%", end="\r")

    def execute(self)->bool:
        try:
            # for choice in self.user_choice:
            logging.info("Download %s microbiome reference catalog", self.choice)
            url = self.catalogues_config[self.choice]["nr_catalogue"]
            md5_expect = self.catalogues_config[self.choice]["md5"]
            catalog_fasta = self.meteor.ref_dir / self.catalogues_config[self.choice]["filename"]
            urlretrieve(url, filename=catalog_fasta, reporthook=self.show_progress)
            assert md5_expect == self.getmd5(catalog_fasta)
            reference_builder = ReferenceBuilder(self.meteor,
                catalog_fasta)
            reference_builder.execute()
        except AssertionError:
            logging.error("MD5sum of %s has a different value than expected", catalog_fasta)
        return True
