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
# from subprocess import run, Popen, PIPE
# from packaging.version import parse

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
from concurrent.futures import ProcessPoolExecutor, as_completed


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

    def set_tree_config(self):
        """Define the census configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "phylogeny": {
                "" "phylogeny_tool": "cogent3",
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

    def process_msp_file(
        self, msp_file: Path, idx: int, msp_count: int, tree_dir, tmp_dir
    ) -> Tuple[Path, bool]:
        """Process a single MSP file and generate its phylogeny tree."""
        logging.info(
            "%d/%d %s: Start analysis",
            idx,
            msp_count,
            msp_file.name.replace(".fasta", ""),
        )
        tree_file = tree_dir / f"{msp_file.stem}.tree"
        # ali_file = tree_dir / f"{msp_file.stem}_aligned.fasta"
        self.tree_files: list[Path] = []
        with NamedTemporaryFile(mode="wt", dir=tmp_dir, suffix=".fasta") as temp_clean:
            # with NamedTemporaryFile(
            #     suffix=".fasta", dir=tmp_dir, delete=True
            # ) as temp_ali_file:
            ## with indel
            # # Start alignment
            # with Popen(
            #     [
            #         "mafft",
            #         "--thread",
            #         str(2),
            #         "--quiet",
            #         str(msp_file.resolve()),
            #     ],
            #     stdout=temp_ali_file,  # Redirect stdout to the temp file
            #     stderr=PIPE,
            # ) as align_exec:
            #     _, error = align_exec.communicate()
            #     if align_exec.returncode != 0:
            #         raise RuntimeError(f"MAFFT failed with error: {error}")
            #     logging.info("Clean sites for %s", msp_file.name)
            #     with ali_file.open("w") as aligned_seq:
            #         _, info_sites = self.clean_sites(
            #             Path(temp_ali_file.name), aligned_seq
            #         )
            # Clean sites
            logging.info("Clean sites")
            _, info_sites = self.clean_sites(msp_file, temp_clean)
            if info_sites < self.min_info_sites:
                logging.info(
                    "Only %d informative sites (< %d threshold) left after cleaning, skipping %s.",
                    info_sites,
                    self.min_info_sites,
                    msp_file.name.replace(".fasta", ""),
                )
                return tree_file, False  # Return False to indicate skipping
            aligned_seqs = load_aligned_seqs(
                temp_clean.name,
                moltype="dna",
            )
            # cleaned_alignment = load_aligned_seqs(ali_file, moltype="dna")
            # d = EstimateDistances(cleaned_alignment, submodel=GTR())
            d = EstimateDistances(aligned_seqs, submodel=GTR())
            d.run(show_progress=False)

            # Create UPGMA Tree
            mycluster = upgma(d.get_pairwise_distances())
            mycluster = mycluster.unrooted_deepcopy()

            with tree_file.open("w") as f:
                f.write(
                    self.remove_edge_labels(mycluster.get_newick(with_distances=True))
                )
        # Perform alignments and UPGMA
        logging.info("Align sequences")

        return tree_file, tree_file.exists()

    def execute(self) -> None:
        logging.info("Launch phylogeny analysis")
        start = perf_counter()
        msp_count = len(self.msp_file_list)
        ## In case of INDEL
        # Check the mafft version
        # mafft_exec = run(["mafft", "--version"], check=False, capture_output=True)
        # if mafft_exec.returncode != 0:
        #     logging.error(
        #         "Checking mafft version failed:\n%s",
        #         mafft_exec.stderr.decode("utf-8"),
        #     )
        #     sys.exit(1)
        # mafft_version = str(mafft_exec.stderr.decode("utf-8")).split(" ")[0][1:]
        # if parse(mafft_version) < self.meteor.MIN_MAFFT_VERSION:
        #     logging.error(
        #         "The mafft version %s is outdated for meteor. Please update mafft to >= %s.",
        #         mafft_version,
        #         self.meteor.MIN_MAFFT_VERSION,
        #     )
        #     sys.exit(1)
        # Using ProcessPoolExecutor to parallelize the MSP file processing
        with ProcessPoolExecutor(max_workers=self.meteor.threads) as executor:
            futures = {
                executor.submit(
                    self.process_msp_file,
                    msp_file,
                    idx,
                    msp_count,
                    self.meteor.tree_dir,
                    self.meteor.tmp_dir,
                ): msp_file
                for idx, msp_file in enumerate(self.msp_file_list, start=1)
            }

            for future in as_completed(futures):
                msp_file = futures[future]
                try:
                    tree_file, success = future.result()
                    if success:
                        self.tree_files.append(tree_file)
                        logging.info(
                            "Completed MSP tree for MSP %s",
                            msp_file.name.replace(".fasta", ""),
                        )
                    else:
                        logging.info(
                            "Skipped MSP %s due to insufficient informative sites",
                            msp_file.name.replace(".fasta", ""),
                        )
                except Exception as exc:
                    logging.error(
                        "MSP %s generated an exception: %s", msp_file.name, exc
                    )

        logging.info("Completed phylogeny in %f seconds", perf_counter() - start)
        logging.info(
            "Trees were generated for %d/%d MSPs", len(self.tree_files), msp_count
        )

        # Save configuration after all trees are processed
        config = self.set_tree_config()
        self.save_config(config, self.meteor.tree_dir / "census_stage_4.json")
