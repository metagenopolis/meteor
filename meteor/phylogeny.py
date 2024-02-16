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

"""Effective phylogeny"""
from dataclasses import dataclass
from meteor.session import Session, Component
from subprocess import call, run
from typing import Type, List, Dict
from time import perf_counter
from pathlib import Path
import tempfile
from tempfile import NamedTemporaryFile
from packaging.version import Version, parse
from collections import OrderedDict
from datetime import datetime
import re
import logging
import sys


@dataclass
class Phylogeny(Session):
    """Run the fastree"""

    meteor: Type[Component]
    msp_file_list: List[Path]
    max_gap: float
    gap_char: str

    def compute_site_info(self, sequences: List[str]) -> List[float]:
        """Calculate the percentage of "_" at each position
        :param sequences: (List) A list of sequence
        :return: (List) A list of the ratio of gap at each position
        """
        return [float(col.count(self.gap_char)) / len(col) for col in zip(*sequences)]

    def clean_sites(
        self, msp_file: Path, output: tempfile._TemporaryFileWrapper
    ) -> None:
        """Clean msp sequence according to a certain level of gap at each sites.
        :param msp_file: (Path) Fasta file
        :param output_file: (Path) Output cleaned fasta file
        """
        gene_dict = OrderedDict(
            (gene_id, seq) for gene_id, seq in self.get_sequences_class(msp_file)
        )
        info_ratio = self.compute_site_info(list(gene_dict.values()))
        for gene_id, seq in gene_dict.items():
            # assert len(info_ratio) == len(seq)
            output_seq = "".join(
                [seq[i] for i, perc in enumerate(info_ratio) if perc <= self.max_gap]
            )
            print(f">{gene_id}\n{output_seq}\n", file=output)

    def set_tree_config(
        self, fasttree_version: str, tree_files: List[Path]
    ) -> Dict:  # pragma: no cover
        """Define the census configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "phylogeny": {
                "phylogeny_tool": "FastTree",
                "phylogeny_version": fasttree_version,
                "phylogeny_date": datetime.now().strftime("%Y-%m-%d"),
                "tree_files": ",".join([tree.name for tree in tree_files]),
            }
        }
        return config

    def execute(self) -> None:
        # Define the regex pattern to match the version number
        version_pattern = re.compile(r"FastTree version (\d+\.\d+\.\d+)")
        fasttree_help = str(run(["FastTree"], capture_output=True).stderr).split("\\n")[
            0
        ]
        # Use findall to extract the version number
        matches = version_pattern.findall(fasttree_help)

        # Check if a match is found
        if matches:
            fasttree_version = matches[0]
            if parse(fasttree_version) < Version("1.9.0"):
                logging.error(
                    "Error, the FastTree version %s is outdated for meteor. Please update FastTree to >=1.9.0.",
                    fasttree_version,
                )
                sys.exit()
        else:
            fasttree_version = "UNK"
            logging.error("Failed to determine the FastTree version.")
        # Start phylogenies
        start = perf_counter()
        tree_files: List[Path] = []
        for msp_file in self.msp_file_list:
            with NamedTemporaryFile(
                mode="wt", dir=self.meteor.tmp_dir, suffix=".fasta", delete=False
            ) as temp_clean:
                tree_file = self.meteor.tree_dir / f"{msp_file.stem}.tree"
                # Clean sites
                self.clean_sites(msp_file, temp_clean)
                with tree_file.open("wt", encoding="UTF-8") as tree:
                    # Compute trees
                    call(
                        [
                            "FastTree",
                            "-nt",
                            "-gtr",
                            temp_clean.name,
                        ],
                        stdout=tree,
                    )
                    tree_files.append(tree_file)
        logging.info("Completed phylogeny in %f seconds", perf_counter() - start)
        config = self.set_tree_config(fasttree_version, tree_files)
        self.save_config(config, self.meteor.tree_dir / f"census_stage_4.json")
