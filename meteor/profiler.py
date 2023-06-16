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

"""Profile the abundance of genes"""

from meteor.session import Session, Component
from dataclasses import dataclass, field
from typing import Type
from configparser import ConfigParser
import pandas as pd

@dataclass
class Profiler(Session):
    """Profile session for abundance and annotation
    """
    meteor: Type[Component]
    normalization: str
    rarefaction_level: int
    statistics: str

    # def normalize_gene_length(self, ):
    #     gene_length = pd.read_csv(, delimiter="\t", header=0)
    #     gene_length.columns = ["gene_id", "gene_name", "Size"]

    # def rarefy(self):
    #     """
    #     """
    #     pd.read_csv(mgs_list, delimiter="\t")

    def execute(self) -> bool:
        """Compute the mapping"""
        # Load gene length
        census_ini_files = list(self.meteor.mapping_dir.glob("**/*_census_stage_1.ini"))
        assert len(census_ini_files) > 0
        # Load counts
        for library in census_ini_files:
            census_ini = ConfigParser()
            with open(library, "rt", encoding="UTF-8") as cens:
                census_ini.read_file(cens)
                # Normalize by gene length
                print(census_ini)
                # Rarefy
                # Combine counts per specie
                # create a matrix

