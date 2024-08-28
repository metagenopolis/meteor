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

# import sys
from dataclasses import dataclass, field
from meteor.session import Session, Component

# from subprocess import check_call, run, DEVNULL
from time import perf_counter
from pathlib import Path
import tempfile
from tempfile import NamedTemporaryFile

# from packaging.version import parse
from collections import OrderedDict
from datetime import datetime
from typing import Iterable, Tuple
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
    min_info_sites: int
    tree_files: list[Path] = field(default_factory=list)

    def compute_site_info(self, sequences: Iterable[str]) -> list[float]:
        """Calculate the percentage of "_" at each position
        :param sequences: (List) A list of sequence
        :return: (List) A list of the ratio of gap at each position
        """
        return [
            float(col.count(self.meteor.DEFAULT_GAP_CHAR)) / len(col)
            for col in zip(*sequences)
        ]

    def clean_sites(
        self, msp_file: Path, output: tempfile._TemporaryFileWrapper
    ) -> Tuple[dict[str, str], int]:
        """Clean msp sequence according to a certain level of gap at each sites.
        :param msp_file: (Path) Fasta file
        :param output_file: (Path) Output cleaned fasta file
        :return: Dict of cleaned sequences
        """
        # Read sequences from msp_file and store in an OrderedDict
        gene_dict = OrderedDict(
            (gene_id, seq) for gene_id, seq in self.get_sequences(msp_file)
        )
        # Compute site information
        info_ratio = self.compute_site_info(gene_dict.values())
        # Count sites with more than the specified maximum gap ratio
        info_sites = sum(1 for ratio in info_ratio if ratio <= self.max_gap)
        logging.info(
            "%d/%d sites with less than %.1f%% gaps",
            info_sites,
            len(info_ratio),
            self.max_gap * 100,
        )
        resultdict = {}
        for gene_id, seq in gene_dict.items():
            # assert len(info_ratio) == len(seq)
            output_seq = "".join(
                seq[i] for i, perc in enumerate(info_ratio) if perc <= self.max_gap
            )
            print(f">{gene_id}\n{output_seq}\n", file=output)
            resultdict[gene_id] = output_seq
        print(flush=True, file=output)
        return resultdict, info_sites
        # return info_sites

    # def set_tree_config(self, raxml_ng_version: str) -> dict:  # pragma: no cover
    def set_tree_config(self):
        """Define the census configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "phylogeny": {
                "phylogeny_tool": "cogent3",
                # "phylogeny_tool": "raxml-ng",
                # "phylogeny_version": raxml_ng_version,
                "phylogeny_date": datetime.now().strftime("%Y-%m-%d"),
                "tree_files": ",".join([tree.name for tree in self.tree_files]),
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
        # raxml_ng_exec = run(["raxml-ng", "--version"], check=False, capture_output=True)
        # if raxml_ng_exec.returncode != 0:
        #     logging.error(
        #         "Checking raxml-ng failed:\n%s", raxml_ng_exec.stderr.decode("utf-8")
        #     )
        #     sys.exit(1)
        # raxml_ng_help = raxml_ng_exec.stdout.decode("utf-8")
        # Define the regex pattern to match the version number
        # version_pattern = re.compile(r"RAxML-NG v\. (\d+\.\d+\.\d+)")
        # match = version_pattern.search(raxml_ng_help)
        # Check if a match is found
        # if not match:
        #     logging.error("Failed to determine the raxml-ng version.")
        #     sys.exit(1)
        # raxml_ng_version = match.group(1)
        # if parse(raxml_ng_version) < self.meteor.MIN_RAXML_NG_VERSION:
        #     logging.error(
        #         "The raxml-ng version %s is outdated for meteor. Please update raxml-ng to >= %s.",
        #         raxml_ng_version,
        #         self.meteor.MIN_RAXML_NG_VERSION,
        #     )
        #     sys.exit(1)

        # Start phylogenies
        start = perf_counter()
        self.tree_files: list[Path] = []
        msp_count = len(self.msp_file_list)
        for idx, msp_file in enumerate(self.msp_file_list, start=1):
            logging.info(
                "%d/%d %s: Start analysis",
                idx,
                msp_count,
                msp_file.name.replace(".fasta", ""),
            )
            with NamedTemporaryFile(
                mode="wt", dir=self.meteor.tmp_dir, suffix=".fasta"
            ) as temp_clean:
                tree_file = self.meteor.tree_dir / f"{msp_file.name}".replace(
                    ".fasta", ""
                )
                # Clean sites
                logging.info("Clean sites")
                _, info_sites = self.clean_sites(msp_file, temp_clean)
                if info_sites < self.min_info_sites:
                    logging.info(
                        "Only %d informative sites (< %d threshold) left after cleaning, skip.",
                        info_sites,
                        self.min_info_sites,
                    )
                # elif len(cleaned_seqs) >= 4:
                #     # Compute trees
                #     logging.info("Run raxml-ng")
                #     result = check_call(
                #         [
                #             "raxml-ng",
                #             "--threads",
                #             "auto{{{}}}".format(self.meteor.threads),
                #             "--workers",
                #             "auto",
                #             "--search1",
                #             "--msa",
                #             temp_clean.name,
                #             "--model",
                #             "GTR+G",
                #             "--redo",
                #             "--force",
                #             "perf,msa",  # not working with raxml-ng-mpi
                #             "--prefix",
                #             str(tree_file.resolve()),
                #         ],
                #         stdout = DEVNULL
                #     )
                # if result != 0:
                #     logging.error("raxml-ng failed with return code %d", result)
                else:
                    # logging.info(
                    #     "Less than 4 sequences, run cogent3"
                    # )
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
                if tree_file.with_suffix(".tree").exists():
                    self.tree_files.append(tree_file.with_suffix(".tree"))
                    logging.info(
                        "Completed MSP tree for MSP %s",
                        msp_file.name.replace(".fasta", ""),
                    )
                # elif tree_file.with_suffix(".raxml.bestTree").exists():
                #     self.tree_files.append(tree_file.with_suffix(".raxml.bestTree"))
                # logging.info("Completed MSP tree with raxml")
                else:
                    logging.info("No tree file generated")
        logging.info("Completed phylogeny in %f seconds", perf_counter() - start)
        logging.info(
            "Trees were generated for %d/%d MSPs", len(self.tree_files), msp_count
        )
        # config = self.set_tree_config(raxml_ng_version)
        config = self.set_tree_config()
        self.save_config(config, self.meteor.tree_dir / "census_stage_4.json")
