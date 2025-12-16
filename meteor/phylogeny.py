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
import pandas as pd

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
from typing import Union


@dataclass
class Phylogeny(Session):
    """Run the fastree"""

    meteor: type[Component]
    msp_file_list: list[Path]
    max_gap: float
    min_info_sites: int
    gtr: bool
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
    ) -> Tuple[dict[Union[int, str, None], str], int]:
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
            "phylogeny": {
                "phylogeny_meteor_version": self.meteor.version,
                "phylogeny_tool": "cogent3",
                "phylogeny_date": datetime.now().strftime("%Y-%m-%d"),
                "tree_mode": "GTR" if self.gtr else "TN93",
                "tree_files": ",".join([tree.name for tree in self.tree_files]),
            },
        }
        return config

    def remove_edge_labels(self, newick: str) -> str:
        # This regular expression matches the edge labels (like "edge.0:")
        pattern = r"\b(edge\.\d+):\b"
        # Replace matched patterns with ":" (effectively removing the edge label)
        return re.sub(pattern, ":", newick)
    
    def _generate_pairwise_comparison_table(self, sequences: dict, dists, output_file: Path) -> None:
        """Generate a pairwise comparison table for each tree with detailed statistics.
        
        Distance categories (based on similarity thresholds):
        - <= 0.0001: same_strain (99.99% similarity)
        - <= 0.01: same_species (≥98% similarity) 
        - <= 0.015: same_subspecies (≥97% similarity)
        - > 0.015: divergent (<97% similarity)
        
        :param sequences: Dictionary of sequence IDs and their sequences or cogent3 Alignment object
        :param dists: Distance matrix from cogent3
        :param output_file: Output TSV file path
        """
        # Define nucleotide information categories
        # Maximal information: strictly A, C, G, or T (unambiguous)
        MAXIMAL_INFO = {'A', 'C', 'G', 'T'}
        
        # Minimal information: A, C, G, T, or IUPAC codes (excluding N, gaps, and ?)
        # IUPAC codes: R, Y, S, W, K, M, B, D, H, V (N is excluded as it's "unknown")
        MINIMAL_INFO = {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}
        
        # Get sequence names and sequence getter function
        if hasattr(sequences, 'keys'):
            seq_names = list(sequences.keys())
            get_seq = lambda name: str(sequences[name])
        else:
            seq_names = list(sequences.names)
            get_seq = lambda name: str(sequences.get_seq(name))
        
        # Pre-convert all sequences to uppercase strings and calculate per-sample statistics
        seq_strs = {}
        per_sample_stats = {}
        for name in seq_names:
            seq_str = get_seq(name).upper()
            seq_strs[name] = seq_str
            total_len = len(seq_str)
            
            # Count positions with minimal and maximal information for this sample
            minimal_info_count = sum(1 for base in seq_str if base in MINIMAL_INFO)
            maximal_info_count = sum(1 for base in seq_str if base in MAXIMAL_INFO)
            
            per_sample_stats[name] = {
                'total_len': total_len,
                'minimal_info_count': minimal_info_count,
                'maximal_info_count': maximal_info_count,
            }
        
        # Prepare data for DataFrame
        table_data = []
        
        # Generate pairwise comparisons (upper triangle only)
        for i, name1 in enumerate(seq_names):
            seq1_str = seq_strs[name1]
            stats1 = per_sample_stats[name1]
            
            for j, name2 in enumerate(seq_names):
                if i >= j:  # Skip diagonal and lower triangle
                    continue
                    
                seq2_str = seq_strs[name2]
                stats2 = per_sample_stats[name2]
                
                # Ensure sequences have the same length
                total_length = stats1['total_len']
                if stats2['total_len'] != total_length:
                    logging.warning(f"Sequence lengths differ for {name1} ({stats1['total_len']}) and {name2} ({stats2['total_len']})")
                    continue
                
                # Initialize overlap counters
                overlap_minimal_info = 0  # Both samples have minimal information (ACGT+IUPAC)
                overlap_maximal_info = 0  # Both samples have maximal information (ACGT only)
                
                # Calculate overlap statistics position by position
                for k in range(total_length):
                    base1 = seq1_str[k]
                    base2 = seq2_str[k]
                    
                    # Check if both positions have minimal information (not N, -, or ?)
                    if base1 in MINIMAL_INFO and base2 in MINIMAL_INFO:
                        overlap_minimal_info += 1
                    
                    # Check if both positions have maximal information (strictly ACGT)
                    if base1 in MAXIMAL_INFO and base2 in MAXIMAL_INFO:
                        overlap_maximal_info += 1
                
                # Calculate percentages
                overlap_minimal_info_pc = (overlap_minimal_info / total_length * 100) if total_length > 0 else 0.0
                overlap_maximal_info_pc = (overlap_maximal_info / total_length * 100) if total_length > 0 else 0.0
                
                minimal_info_pc_sample1 = (stats1['minimal_info_count'] / total_length * 100) if total_length > 0 else 0.0
                maximal_info_pc_sample1 = (stats1['maximal_info_count'] / total_length * 100) if total_length > 0 else 0.0
                
                minimal_info_pc_sample2 = (stats2['minimal_info_count'] / total_length * 100) if total_length > 0 else 0.0
                maximal_info_pc_sample2 = (stats2['maximal_info_count'] / total_length * 100) if total_length > 0 else 0.0
                
                # Get distance from distance matrix
                if hasattr(dists, 'get_distance'):
                    distance = dists.get_distance(name1, name2)
                else:
                    distance = dists[i, j] if hasattr(dists, '__getitem__') else 0.0
                
                # Categorize distance based on similarity thresholds
                if distance <= 0.0001:
                    distance_category = "same_strain"
                elif distance <= 0.01:
                    distance_category = "same_species"
                elif distance <= 0.015:
                    distance_category = "same_subspecies"
                else:
                    distance_category = "divergent"
                
                # Add row to table data with the requested column names
                table_data.append({
                    'sample1': name1,
                    'sample2': name2,
                    'total_length': total_length,
                    'overlap_noN_info_count': overlap_minimal_info,
                    'overlap_noIUPAC_info_count': overlap_maximal_info,
                    'overlap_noN_info_pc': overlap_minimal_info_pc,
                    'overlap_noIUPAC_info_pc': overlap_maximal_info_pc,
                    'noN_info_pc_sample1': minimal_info_pc_sample1,
                    'noN_info_pc_sample2': minimal_info_pc_sample2,
                    'noIUPAC_info_pc_sample1': maximal_info_pc_sample1,
                    'noIUPAC_info_pc_sample2': maximal_info_pc_sample2,
                    'distance': distance,
                    'distance_category': distance_category,
                })
        
        # Create DataFrame and save to TSV
        if table_data:
            df = pd.DataFrame(table_data)
            df.to_csv(output_file, sep='\t', index=False)
            logging.info(f"Pairwise comparison table saved to {output_file}")
        else:
            # Create empty file with headers if no data
            columns = ['sample1', 'sample2', 'total_length', 'overlap_noN_info_count',
                    'overlap_noIUPAC_info_count', 'overlap_noN_info_pc', 'overlap_noIUPAC_info_pc',
                    'noN_info_pc_sample1', 'noN_info_pc_sample2', 'noIUPAC_info_pc_sample1',
                    'noIUPAC_info_pc_sample2', 'distance', 'distance_category']
            pd.DataFrame(columns=columns).to_csv(output_file, sep='\t', index=False)
            logging.info(f"Empty pairwise comparison table saved to {output_file}")
    # def _write_distance_matrix_to_tsv(self, dists, dist_file: Path):
    #     """Write distance matrix to TSV file."""        
    #     # Get sequence names
    #     seq_names = list(dists.names)
        
    #     # Convert distance matrix to DataFrame
    #     dist_data = []
    #     for i, name1 in enumerate(seq_names):
    #         row = []
    #         for j, name2 in enumerate(seq_names):
    #             if hasattr(dists, 'get_distance'):
    #                 # For cogent3 distance matrices
    #                 distance = dists.get_distance(name1, name2)
    #             else:
    #                 # For numpy-like matrices
    #                 distance = dists[i, j]
    #             row.append(distance)
    #         dist_data.append(row)
        
    #     # Create DataFrame and save to TSV
    #     df = pd.DataFrame(dist_data, index=seq_names, columns=seq_names)
    #     df.to_csv(dist_file, sep='\t')
        
    #     logging.info(f"Distance matrix saved to {dist_file}")

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
        dist_file = tree_dir / f"{msp_file.stem}.tsv"
        tree_file = tree_dir / f"{msp_file.stem}.tree"
        # Check if tree file already exists and is not empty
        if tree_file.exists() and tree_file.stat().st_size > 0:
            # Also check if distance file exists (optional but recommended)
            if dist_file.exists() and dist_file.stat().st_size > 0:
                logging.info(
                    "%d/%d %s: Tree file already exists, skipping computation",
                    idx,
                    msp_count,
                    msp_file.name.replace(".fasta", ""),
                )
                return tree_file, True
            else:
                logging.warning(
                    "%d/%d %s: Tree file exists but distance file is missing, reprocessing",
                    idx,
                    msp_count,
                    msp_file.name.replace(".fasta", ""),
                )
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
                format_name="fasta",
            )
            # cleaned_alignment = load_aligned_seqs(ali_file, moltype="dna")
            # d = EstimateDistances(cleaned_alignment, submodel=GTR())
            if self.gtr:
                d = EstimateDistances(aligned_seqs, submodel=GTR())
                d.run(show_progress=False)
                # Create UPGMA Tree
                dists = d.get_pairwise_distances()
                mycluster = upgma(dists)
            else:
                dists = aligned_seqs.distance_matrix(calc="tn93")
                mycluster = upgma(dists)

            # Save distance matrix
            # self._write_distance_matrix_to_tsv(dists, dist_file)

            mycluster = mycluster.unrooted_deepcopy()
            # Get distance matrix and tip order
            ultrametric_dist = mycluster.tip_to_tip_distances()
 
            # Convert to squareform DataFrame
            df = pd.DataFrame(ultrametric_dist, index=ultrametric_dist.names, columns=ultrametric_dist.names)
            df.to_csv(dist_file, sep="\t")

            # Generate pairwise comparison table
            comparison_table_file = tree_dir / f"{msp_file.stem}_comparison.tsv"
            self._generate_pairwise_comparison_table(aligned_seqs, dists, comparison_table_file)

            with tree_file.open("w") as f:
                f.write(
                    self.remove_edge_labels(mycluster.get_newick(with_distances=True))
                )

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
