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
import re
import logging
import sys
from dataclasses import dataclass
from meteor.session import Session, Component
from subprocess import check_call, run
from time import perf_counter
from pathlib import Path
import tempfile
from tempfile import NamedTemporaryFile
from packaging.version import parse
from collections import OrderedDict
from datetime import datetime
from typing import Iterable
from cogent3 import load_aligned_seqs
from cogent3.evolve.distance import EstimateDistances
from cogent3.evolve.models import GTR
from cogent3.cluster.UPGMA import upgma


@dataclass
class Phylogeny(Session):
    """Run the fastree"""

    meteor: type[Component]
    msp_file_list: list[Path]
    max_gap: float
    gap_char: str

    def compute_site_info(self, sequences: Iterable[str]) -> list[float]:
        """Calculate the percentage of "_" at each position
        :param sequences: (List) A list of sequence
        :return: (List) A list of the ratio of gap at each position
        """
        return [float(col.count(self.gap_char)) / len(col) for col in zip(*sequences)]

    def clean_sites(
        self, msp_file: Path, output: tempfile._TemporaryFileWrapper
    ) -> dict[str, str]:
        """Clean msp sequence according to a certain level of gap at each sites.
        :param msp_file: (Path) Fasta file
        :param output_file: (Path) Output cleaned fasta file
        :return: Dict of cleaned sequences
        """
        gene_dict = OrderedDict(
            (gene_id, seq) for gene_id, seq in self.get_sequences_class(msp_file)
        )
        info_ratio = self.compute_site_info(gene_dict.values())
        resultdict = {}
        for gene_id, seq in gene_dict.items():
            # assert len(info_ratio) == len(seq)
            output_seq = "".join(
                seq[i] for i, perc in enumerate(info_ratio) if perc <= self.max_gap
            )
            print(f">{gene_id}\n{output_seq}\n", file=output)
            resultdict[gene_id] = output_seq
        return resultdict

    def compute_mutation_rate(self, seq1, seq2):
        """Compute mutation rate between two sequences."""
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length.")

        mutations = sum(nuc1 != nuc2 for nuc1, nuc2 in zip(seq1, seq2))
        return mutations / len(seq1)

    def set_tree_config(
        self, raxml_ng_version: str, tree_files: list[Path]
    ) -> dict:  # pragma: no cover
        """Define the census configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "phylogeny": {
                "phylogeny_tool": "raxml-ng",
                "phylogeny_version": raxml_ng_version,
                "phylogeny_date": datetime.now().strftime("%Y-%m-%d"),
                "tree_files": ",".join([tree.name for tree in tree_files]),
            },
        }
        return config

    def remove_edge_labels(self, newick: str) -> str:
        # This regular expression matches the edge labels (like "edge.0:")
        pattern = r"\b(edge\.\d+):\b"
        # Replace matched patterns with ":" (effectively removing the edge label)
        return re.sub(pattern, ":", newick)

    def execute(self) -> None:
        logging.info("Launch phylogeny analysis")
        raxml_ng_exec = run(["raxml-ng", "--version"], check=False, capture_output=True)
        if raxml_ng_exec.returncode != 0:
            logging.error(
                "Checking raxml-ng failed:\n%s", raxml_ng_exec.stderr.decode("utf-8")
            )
            sys.exit(1)
        raxml_ng_help = raxml_ng_exec.stdout.decode("utf-8")
        # Define the regex pattern to match the version number
        version_pattern = re.compile(r"RAxML-NG v\. (\d+\.\d+\.\d+)")
        match = version_pattern.search(raxml_ng_help)
        # Check if a match is found
        if not match:
            logging.error("Failed to determine the raxml-ng version.")
            sys.exit(1)
        raxml_ng_version = match.group(1)
        if parse(raxml_ng_version) < self.meteor.MIN_RAXML_NG_VERSION:
            logging.error(
                "The raxml-ng version %s is outdated for meteor. Please update raxml-ng to >= %s.",
                raxml_ng_version, self.meteor.MIN_RAXML_NG_VERSION
            )
            sys.exit(1)
            
        # Start phylogenies
        start = perf_counter()
        tree_files: list[Path] = []
        msp_count = len(self.msp_file_list)
        for idx, msp_file in enumerate(self.msp_file_list, start=1):
            logging.info(
                "Start analysis of MSP %s: %d/%d", msp_file.name, idx, msp_count
            )
            with NamedTemporaryFile(
                mode="wt", dir=self.meteor.tmp_dir, suffix=".fasta"
            ) as temp_clean:
                tree_file = self.meteor.tree_dir / f"{msp_file.name}".replace(
                    ".fasta", ""
                )
                # Clean sites
                cleaned_seqs = self.clean_sites(msp_file, temp_clean)
                logging.info("Clean sites for MSP %d/%d", idx, msp_count)
                if len(cleaned_seqs) >= 4:
                    # Compute trees
                    result = check_call(
                        [
                            "raxml-ng",
                            "--threads",
                            str(self.meteor.threads),
                            "--msa",
                            temp_clean.name,
                            "--model",
                            "GTR+G",
                            "--redo",
                            # "--force perf_threads", # only with raxml-ng-mpi
                            "--prefix",
                            str(tree_file.resolve()),
                        ]
                    )
                    if result != 0:
                        logging.error("raxml-ng failed with return code %d", result)
                else:
                    logging.info(
                        "MSP %s have less than 4 sequences, distance will be calculated with cogent3",
                        msp_file.name,
                    )
                    aligned_seqs = load_aligned_seqs(
                        temp_clean.name,
                        moltype="dna",
                    )
                    d = EstimateDistances(aligned_seqs, submodel=GTR())
                    d.run(show_progress=False)
                    mycluster = upgma(d.get_pairwise_distances())
                    mycluster = mycluster.unrooted_deepcopy()
                    with tree_file.with_suffix(".tree").open("w") as f:
                        f.write(
                            self.remove_edge_labels(
                                mycluster.get_newick(with_distances=True)
                            )
                        )
                    # Edges get a name which is not supported by ete3
                    # mycluster.write(
                    #     tree_file.with_suffix(".tree"),
                    # )
                tree_files.append(tree_file)
            logging.info("Completed MSP tree %d/%d", idx, msp_count)
        logging.info("Completed phylogeny in %f seconds", perf_counter() - start)
        config = self.set_tree_config(raxml_ng_version, tree_files)
        self.save_config(config, self.meteor.tree_dir / "census_stage_4.json")
