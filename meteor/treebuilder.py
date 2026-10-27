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

"""Process msp and compute tree"""
import logging
import sys
import lzma
# import ete3  # type: ignore[import]
import os
import fnmatch
from pathlib import Path
from collections import defaultdict
from meteor.session import Session, Component
from meteor.phylogeny import Phylogeny
from dataclasses import dataclass
from tempfile import mkdtemp
from ete4 import Tree
from shutil import rmtree
from concurrent.futures import ThreadPoolExecutor
from typing import ClassVar


@dataclass
class TreeBuilder(Session):
    """Counter session map and count"""

    DEFAULT_MAX_GAP: ClassVar[float] = 1.0
    OUTPUT_FORMATS: ClassVar[list[str | None]] = [None, "png", "svg", "pdf", "txt"]
    DEFAULT_OUTPUT_FORMAT: ClassVar[str | None] = None
    DEFAULT_WIDTH: ClassVar[int] = 500
    DEFAULT_HEIGHT: ClassVar[int] = 500
    DEFAULT_NUM_THREADS: ClassVar[int] = 1
    DEFAULT_MIN_INFO_SITES: ClassVar[int] = 10

    meteor: type[Component]
    max_gap: float
    min_info_sites: int
    gtr: bool
    width: int
    height: int
    format: str | None
    tshrink: bool = False
    sample_file: Path | None = None
    # Outputs of an MSP tree, removed when its samples change
    TREE_OUTPUTS: ClassVar[list[str]] = [
        ".tree", ".tsv", "_comparison.tsv", "_qcflags.tsv", "_shrink.tsv",
    ]

    def __post_init__(self) -> None:
        if self.format not in TreeBuilder.OUTPUT_FORMATS:
            raise ValueError(f"{self.format} is not a valid output format")

        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.tree_dir.mkdir(exist_ok=True, parents=True)

    # def concatenate(self, msp_file_dict: dict[str, list[Path]]) -> list[Path]:
    #     """Concatenate fasta file from distinct samples
    #     :param msp_file_dict: (dict) A dictionnary of each msp version in the different samples
    #     :return: (list) A list of all concatenated fasta
    #     """
    #     msp_list = []
    #     # Concatenate files that occur in more than one directory
    #     for filename, paths in msp_file_dict.items():
    #         if len(paths) > 1:
    #             res = self.meteor.tree_dir / f"{filename}".replace(
    #                 ".fasta.xz", ".fasta"
    #             )
    #             with res.open("wt", encoding="UTF-8") as outfile:
    #                 for path in paths:
    #                     with lzma.open(path, "rt") as infile:
    #                         outfile.write(infile.read())
    #             msp_list += [res]
    #     logging.info("%d MSPs are available for tree analysis.", len(msp_list))
    #     return msp_list

    def concatenate(
        self,
        msp_file_dict: dict[str, list[Path]],
        sample_names: dict[Path, str] | None = None,
    ) -> list[Path]:
        """Concatenate fasta file from distinct samples
        :param msp_file_dict: (dict) A dictionary of each msp version in the different samples
        :param sample_names: {sample directory: sample name} (default: directory name)
        :return: (list) A list of all concatenated fasta
        """
        sample_names = sample_names or {}
        msp_list = []
        to_concatenate = []
        for filename, paths in msp_file_dict.items():
            if len(paths) > 1:
                res = self.meteor.tree_dir / f"{filename}".replace(".fasta.xz", ".fasta")
                
                # Reuse the file of a previous run if it has the same samples
                if res.exists() and res.stat().st_size > 0:
                    logging.debug("File %s exists, checking if concatenation is needed...", res)
                    expected = {sample_names.get(path.parent, path.parent.name) for path in paths}
                    if self.fasta_samples(res) == expected:
                        msp_list.append(res)
                        continue
                    # samples changed (sample list, new samples): recompute the tree
                    logging.info("%s: samples changed, the tree is recomputed.", res.stem)
                    for suffix in self.TREE_OUTPUTS:
                        (self.meteor.tree_dir / f"{res.stem}{suffix}").unlink(missing_ok=True)
                to_concatenate.append((res, paths))
                msp_list.append(res)

        def concatenate_one(res: Path, paths: list[Path]) -> None:
            # Concatenate if file doesn't exist or is empty
            logging.debug("Creating concatenated file %s", res)
            with res.open("wt", encoding="UTF-8") as outfile:
                for path in paths:
                    with lzma.open(path, "rt") as infile:
                        outfile.write(infile.read())

        # Files are independent: decompress them in parallel (lzma releases the GIL)
        with ThreadPoolExecutor(max_workers=max(1, self.meteor.threads)) as executor:
            for future in [executor.submit(concatenate_one, res, paths) for res, paths in to_concatenate]:
                future.result()

        logging.info("%d MSPs are available for tree analysis.", len(msp_list))
        return msp_list

    # def get_msp_distance(self, tree: ete3.TreeNode) -> pd.DataFrame:
    #     samples = list(tree)
    #     distance_matrix = pd.DataFrame(
    #         index=[n.name for n in samples], columns=[n.name for n in samples]
    #     )
    #     distance_matrix.values[[range(len(samples))] * 2] = 0.0
    #     for i in samples[:-1]:
    #         for j in samples[1:]:
    #             distance = tree.get_distance(i, j)
    #             distance_matrix.loc[i.name, j.name] = distance
    #             distance_matrix.loc[j.name, i.name] = distance
    #     return distance_matrix

    @staticmethod
    def fasta_samples(fasta_file: Path) -> set[str]:
        """Sample names (sequence names) of a concatenated MSP fasta"""
        with fasta_file.open("rt", encoding="UTF-8") as fasta:
            return {line[1:].strip() for line in fasta if line.startswith(">")}

    def execute(self) -> None:
        "Merge all files generated by either mapper or profiler."
        # Fetch all census ini files
        #all_census = list(Path(self.meteor.strain_dir).glob("**/*census_stage_3.json"))
        all_census = [Path(root) / f for root, _, files in os.walk(self.meteor.strain_dir, followlinks=True)
                      for f in fnmatch.filter(files, "*census_stage_3.json")]
        if len(all_census) == 0:
            logging.error(
                "No *_census_stage_3.json files found in %s",
                self.meteor.strain_dir,
            )
            sys.exit(1)
        else:
            logging.info("%d samples have been detected.", len(all_census))
        # {sample name: sample directory}, restricted to the sample list (-l)
        sample_dirs = {
            self.read_json(census)["sample_info"]["sample_name"]: census.parent
            for census in all_census
        }
        selected_dirs = None
        if self.sample_file is not None:
            selected_dirs = set(
                self.select_samples(sample_dirs, self.sample_file).values()
            )
        msp_file_dict = defaultdict(list)
        # for filepath in self.meteor.strain_dir.glob("**/*.fasta.xz"):
            # if not filepath.name.endswith("_consensus.fasta.xz"):
                # msp_file_dict[filepath.name].append(filepath)
        for root, _, files in os.walk(self.meteor.strain_dir, followlinks=True):
            if selected_dirs is not None and Path(root) not in selected_dirs:
                continue
            for filename in fnmatch.filter(files, "*.fasta.xz"):
                if not filename.endswith("_consensus.fasta.xz"):
                    msp_file_dict[filename].append(Path(root) / filename)

        # Concatenate msp files
        msp_file_list = self.concatenate(
            msp_file_dict, {path: name for name, path in sample_dirs.items()}
        )
        # Compute phylogenies
        phylogeny_process = Phylogeny(
            self.meteor,
            msp_file_list,
            self.max_gap,
            self.min_info_sites,
            self.gtr,
            tshrink=self.tshrink,
        )
        phylogeny_process.execute()
        # Analyze tree data
        for tree_file in phylogeny_process.tree_files:
            output = tree_file
            if tree_file.suffix == ".bestTree":
                output = tree_file.with_suffix("")
            img_file = self.meteor.tree_dir / f"{output.stem}.{self.format}"
            # try:
            # # Generate a distance msp by msp
            # matrix = self.get_msp_distance(msp_tree)
            # matrix.to_csv(
            #     self.meteor.tree_dir / f"{output.stem}.tsv",
            #     sep="\t",
            #     quoting=csv.QUOTE_NONE,
            # )
            # Draw trees
            # Only load tree if we actually need to use it
            if self.format:  # Only proceed if format is specified
                msp_tree = Tree(str(tree_file.resolve()))
            
                if self.format == "txt":
                    with img_file.open("wt", encoding="UTF-8") as outfile:
                        outfile.write(msp_tree.get_ascii(show_internal=True))
                else:
                    # ts = TreeStyle()
                    # ts.show_leaf_name = True
                    # ts.show_branch_length = True
                    msp_tree.render(
                        str(img_file.resolve()),
                        w=self.width,
                        h=self.height,
                        units="px",
                        dpi=300,
                    )
            # except ete3.parser.newick.NewickError:
            #     logging.info("Not sufficient info in %s.", str(tree_file.resolve()))
        rmtree(self.meteor.tmp_dir, ignore_errors=True)
