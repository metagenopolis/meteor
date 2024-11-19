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
from meteor.parser import Parser
from dataclasses import dataclass, field
import pandas as pd
from pathlib import Path
import importlib.resources
import numpy as np
import logging
import sys
import lzma
from datetime import datetime
from typing import ClassVar


@dataclass
class Profiler(Session):
    """Profile session for abundance and annotation"""

    NO_RAREFACTION: ClassVar[int] = 0
    DEFAULT_RAREFACTION_LEVEL: ClassVar[int] = NO_RAREFACTION
    DEFAULT_RANDOM_SEED: ClassVar[int] = 1234
    NORMALIZATIONS: ClassVar[list[str]] = ["coverage", "fpkm", "raw"]
    DEFAULT_NORMALIZATION: ClassVar[str] = "coverage"
    DEFAULT_COVERAGE_FACTOR: ClassVar[float] = 100.0
    DEFAULT_MSP_FILTER_COMPLETE: ClassVar[float] = 0.1
    DEFAULT_MSP_FILTER_TAXO: ClassVar[float] = 0.2
    DEFAULT_COMPLETENESS: ClassVar[float] = 0.9

    meteor: type[Component]
    rarefaction_level: int
    seed: int
    normalization: str | None
    core_size: int
    msp_filter_user: float | None
    completeness: float
    coverage_factor: float
    msp_filter: float = field(default_factory=float)

    def __post_init__(self):
        if self.normalization not in Profiler.NORMALIZATIONS:
            raise ValueError(f"{self.normalization} is not a valid normalization")

        # Get the json file
        self.sample_config = self.get_census_stage(self.meteor.mapping_dir, 1)

        # Add session info
        config_session = {}
        config_session.update(
            {
                "meteor_version": self.meteor.version,
                "profiling_date": str(datetime.now().strftime("%Y-%m-%d")),
            }
        )
        self.sample_config = self.update_json(
            self.sample_config, "profiling_session", config_session
        )

        # Initialize the unmapped_count
        try:
            self.total_reads = int(self.sample_config["mapping"]["total_read_count"])
            self.mapped_reads = int(self.sample_config["counting"]["counted_reads"])
            self.unmapped_reads = self.total_reads - self.mapped_reads
            logging.info("Number of unmapped reads: %s", str(self.unmapped_reads))
        except KeyError:
            logging.info("No unmapped reads information found. Set to 0.")
            self.unmapped_reads = 0

        # Initialize the sample name
        self.sample_name = self.sample_config["sample_info"]["sample_name"]

        # Prepare the census_stage_2 directory
        self.stage2_dir = self.meteor.profile_dir / self.sample_name
        self.stage2_dir.mkdir(exist_ok=True, parents=True)

        # Initialize the ini ref config parser:
        self.ref_config = self.get_reference_info(self.meteor.ref_dir)

        # Get the database type
        self.database_type = self.ref_config["reference_info"]["database_type"]

        if not self.msp_filter_user:
            if self.database_type == "complete":
                self.msp_filter = Profiler.DEFAULT_MSP_FILTER_COMPLETE
            else:
                self.msp_filter = Profiler.DEFAULT_MSP_FILTER_TAXO
        else:
            self.msp_filter = self.msp_filter_user

        # Get the associated count table
        self.input_count_table = self.meteor.mapping_dir / f"{self.sample_name}.tsv.xz"
        try:
            assert self.input_count_table.is_file()
        except AssertionError:
            logging.error("The count table %s does not exist.", self.input_count_table)
            sys.exit(1)

        # Add a symlink to get the raw count table in the profile directory (for merging purpose)
        raw_count_table_symlink = self.stage2_dir / f"{self.sample_name}_raw.tsv.xz"
        try:
            raw_count_table_symlink.symlink_to(self.input_count_table.resolve())
        except FileExistsError:
            logging.info(
                "The symbolic link to raw data already exists in the profile directory."
            )

        # Check the input count table
        # self.check_file(
        #     self.input_count_table,
        #     {
        #         "gene_id",
        #         "value",
        #         "gene_length",
        #     },
        # )

        # Load the count table
        self.gene_count = self.load_data(self.input_count_table)
        self.gene_count["value"] = self.gene_count["value"].astype(
            pd.SparseDtype("float", fill_value=0.0)
        )

        # Define output filenames
        self.output_base_filename = f"{self.sample_name}"

        # Define MSP filename
        self.msp_filename = (
            self.meteor.ref_dir
            / self.ref_config["reference_file"]["database_dir"]
            / self.ref_config["annotation"]["msp"]["filename"]
        )
        assert self.msp_filename.is_file()
        # self.check_file(
        #     self.msp_filename,
        #     {
        #         "msp_name",
        #         "gene_id",
        #         "gene_category",
        #     },
        # )

        # Get functional db filenames
        if self.database_type == "complete":
            self.db_filenames = {}
            for db in [
                "kegg",
                "mustard",
                "dbcan",
                "eggnog",
                "tigrfam",
                "resfinder",
                "resfinderfg",
            ]:
                self.db_filenames[db] = (
                    self.meteor.ref_dir
                    / self.ref_config["reference_file"]["database_dir"]
                    / self.ref_config["annotation"][db]["filename"]
                )
                assert self.db_filenames[db].is_file()

            # Initialize the module definition file
            self.module_path = (
                importlib.resources.files("meteor") / "data/modules_definition.feather"
            )
            assert self.module_path.is_file()

    def rarefy(self, rarefaction_level: int, unmapped_reads: int, seed: int) -> None:
        """Perform rarefaction on the gene count table.

        :param rarefaction_level: number of reads that will be randomly selected.
        :param unmapped_reads: number of unmapped_reads that should be taken into account
        :param seed: seed to reproduce the random selection of reads.
        """
        try:
            assert rarefaction_level > Profiler.NO_RAREFACTION
        except AssertionError:
            logging.error("You are trying to rarefy with a null or negative number.")
            sys.exit(1)
        # Rarefaction must be performed on integer values
        self.gene_count["value"] = self.gene_count["value"].astype(float).round(0)
        self.gene_count["value"] = self.gene_count["value"].astype(int)
        # Add the unmapped count
        self.gene_count.loc[len(self.gene_count)] = {
            "gene_id": -1,
            "gene_length": 1000,
            "value": unmapped_reads,
        }

        # Check if rarefaction is possible
        if self.gene_count["value"].sum() > rarefaction_level:
            # Transform into a long array of gene_id according to gene occurence
            array_to_rarefy = np.repeat(
                self.gene_count["gene_id"], self.gene_count["value"]
            )
            # Randomly choose among the long array the selected reads
            rng = np.random.default_rng(seed=seed)
            array_rarefied = rng.choice(
                array_to_rarefy, size=rarefaction_level, replace=False
            )
            # Count gene_id occurence to get back to short gene_id list
            unique, counts = np.unique(array_rarefied, return_counts=True)
            self.gene_count["value"] = np.zeros(self.gene_count.shape[0], dtype="int")
            self.gene_count.loc[
                self.gene_count.reset_index()
                .set_index("gene_id")
                .loc[unique, "index"]
                .values,
                "value",
            ] = counts

        # Remove the counts for the gene "-1" (unmapped_reads)
        self.gene_count = self.gene_count[self.gene_count["gene_id"] != -1]

    def normalize_coverage(self, trim_length: int) -> None:
        """Normalize gene count table by coverage.

        :param trim_length: length used in bowtie2 to trim reads before mapping.
        """
        condlist = [
            self.gene_count["gene_length"] >= 2 * trim_length,
            self.gene_count["gene_length"] % 2 == 0,
            self.gene_count["gene_length"] % 2 == 1,
        ]

        choicelist = [
            (self.coverage_factor * self.gene_count["value"])
            / (self.gene_count["gene_length"] - trim_length + 1),
            (self.coverage_factor * self.gene_count["value"] * 4 * trim_length)
            / (self.gene_count["gene_length"] * (self.gene_count["gene_length"] + 2)),
            (self.coverage_factor * self.gene_count["value"] * 4 * trim_length)
            / (self.gene_count["gene_length"] + 1) ** 2,
        ]
        self.gene_count["value"] = np.select(condlist, choicelist)

    def normalize_fpkm(self, rarefaction_level: int, unmapped_reads: int) -> None:
        """Normalize matrix using fpkm method

        :param rarefaction_level: Value of rarefaction level
        :param unmapped_reads: Value of the unmapped reads
        """
        # Compute the unmapped reads after rarefaction
        if rarefaction_level > Profiler.NO_RAREFACTION:
            # Force to 0 if the result is negative (ie, no rarefaction was performed)
            unmapped_reads_after_rf = min(
                rarefaction_level - self.gene_count["value"].sum(), unmapped_reads
            )
        else:
            unmapped_reads_after_rf = unmapped_reads
        # Add the unmapped reads after rf to the gene count table as a pseudo gene
        self.gene_count.loc[len(self.gene_count)] = {
            "gene_id": -1,
            "gene_length": 1000,
            "value": unmapped_reads_after_rf,
        }
        # Normalize
        self.gene_count["value"] = (
            self.gene_count["value"] / self.gene_count["gene_length"]
        )
        self.gene_count["value"] = (
            self.gene_count["value"] / self.gene_count["value"].sum()
        )

        # Remove the counts for the gene "-1" (unmapped_reads)
        self.gene_count = self.gene_count[self.gene_count["gene_id"] != -1]

    def compute_msp(self, msp_dict: dict[str, set[str]], filter_pc: float) -> None:
        """Compute msp abundance table.

        :param msp_dict: dictionnary of msp definition {'msp1': {'gene1', 'gene2'}}
        :param filter_pc: minimum ratio of msp genes that should be detected, otherwise
        msp abundance is set to 0.
        """
        # Compute how many genes are seen for a given msp
        # Restrict to gene table to core genes
        all_core_genes = {item for sublist in msp_dict.values() for item in sublist}
        gene_count_core = self.gene_count.loc[
            self.gene_count["gene_id"].isin(all_core_genes)
        ]
        msp_filter = {
            msp: (
                gene_count_core.loc[
                    gene_count_core["gene_id"].isin(set_genes),
                    "value",
                ]
                > 0
            ).sum()
            / len(set_genes)
            for (msp, set_genes) in msp_dict.items()
        }
        # Compute mean abundance if gene count is above filter threshold, otherwise 0
        msp_table_dict = {
            msp: (
                gene_count_core.loc[
                    gene_count_core["gene_id"].isin(set_genes),
                    "value",
                ].mean()
                if msp_filter[msp] >= filter_pc
                else 0
            )
            for (msp, set_genes) in msp_dict.items()
        }
        self.msp_table = (
            pd.DataFrame.from_dict(msp_table_dict, orient="index", columns=["value"])
            .reset_index()
            .rename(columns={"index": "msp_name"})
        )

    def get_msp_core(self, msp_def_filename: Path, core_size: int) -> dict:
        """Get genes for each msp that correspond to gene core, in the limit of the core_size.

        :param msp_def_filename: path to the msp definition file
        :param core_size: maximum number of core genes to consider.
        """
        # Load msp file
        msp_df = self.load_data(msp_def_filename)
        # Restrict to core
        msp_df_selection = msp_df.loc[msp_df["gene_category"] == "core"]
        # Return the df as a dict of set
        msp_dict = (
            msp_df_selection.groupby("msp_name")["gene_id"]
            .apply(lambda x: set(x.head(core_size)))
            .to_dict()
        )
        return msp_dict

    def compute_msp_stats(self, msp_def_filename: Path) -> float:
        """Compute percentage of reads that map on genes from MSP.

        :param msp_def_filename: A path object pointing to an MSP definition file.
        """
        # Load msp file
        msp_df = self.load_data(msp_def_filename)
        # Get the ensemble of genes used in MSP
        all_msp_genes = msp_df["gene_id"].unique()
        # Get the percentage of reads that map on an MSP
        msp_reads_pc = (
            self.gene_count.loc[
                self.gene_count["gene_id"].isin(all_msp_genes),
                "value",
            ].sum()
            / self.gene_count["value"].sum()
        )
        return round(msp_reads_pc, 2)

    def compute_ko_abundance(self, annot_file: Path) -> None:
        """Compute abundance of enzyme as sum of genes assigned to this enzyme.

        :param annot_file: a path object pointing to the annotation gene_name -> enzyme file.
        """
        # Load annotation file
        annot_df = self.load_data(annot_file)
        # Merge count table and gene annotation
        merged_df = pd.merge(
            annot_df,
            self.gene_count,
            left_on="gene_id",
            right_on="gene_id",
        )
        # Compute sum of KO
        aggregated_count = merged_df.groupby("annotation")["value"].sum().reset_index()
        self.functions = aggregated_count

    def compute_ko_abundance_by_msp(
        self, annot_file: Path, msp_def_filename: Path
    ) -> None:
        """Compute abundance of enzyme as the sum of MSP carrying the enzyme.

        :param annot_file: a path object pointing to the annotation gene_name -> enzyme file.
        :param msp_def_filename: A path object pointing to an MSP definition file.
        """
        # Load annotation file
        annot_df = self.load_data(annot_file)
        # Get KO list
        all_ko = annot_df["annotation"].unique()
        # Load MSP file
        msp_df = self.load_data(msp_def_filename)
        # Merge both data frames
        msp_df_annotated = pd.merge(msp_df, annot_df)
        # Restrict to detected genes
        detected_genes = self.gene_count.loc[self.gene_count["value"] > 0, "gene_id"]
        msp_df_annotated = msp_df_annotated.loc[
            msp_df_annotated["gene_id"].isin(detected_genes)
        ]
        # Create a dict ko: {msp1, msp2}
        ko_dict = (
            msp_df_annotated.groupby("annotation")["msp_name"].apply(set).to_dict()
        )
        # Loop on all_ko since some ko have no msp or ne detected genes
        ko_dict_ab = {
            ko: (
                self.msp_table.loc[
                    self.msp_table["msp_name"].isin(ko_dict[ko]),
                    "value",
                ].sum()
                if ko in ko_dict
                else 0
            )
            for ko in all_ko
        }
        self.functions = (
            pd.DataFrame.from_dict(ko_dict_ab, orient="index", columns=["value"])
            .reset_index()
            .rename(columns={"index": "annotation"})
        )

    def compute_ko_stats(
        self, annot_file: Path, by_msp: bool, msp_def_filename: Path
    ) -> float:
        """Compute percentage of reads that map on annotated genes.

        :param annot_file: path to the annotation file (kegg, mustard, etc)
        :param by_msp: should the genes be restricted to those belonging to an MSP
        :param msp_def_filename: A path object pointing to an MSP definition file.
        """
        # Load annotation file
        annot_df = self.load_data(annot_file)
        if by_msp:
            # Load MSP file
            msp_df = self.load_data(msp_def_filename)
            # Merge both data frames
            annot_df = pd.merge(msp_df, annot_df)
        # Get the genes in MSP AND annotated
        all_msp_genes = annot_df["gene_id"].unique()
        # Get the percentage of reads that map on these genes
        annot_reads_pc = (
            self.gene_count.loc[
                self.gene_count["gene_id"].isin(all_msp_genes),
                "value",
            ].sum()
            / self.gene_count["value"].sum()
        )
        return round(annot_reads_pc, 2)

    def merge_catalogue_info(
        self, msp_file: Path, annot_file: dict[str, Path]
    ) -> pd.DataFrame:
        """Merge info from msp definition and genes functional annotation.

        :param msp_file: path to the msp definition file
        :param annot_file: path to the gene functional annotation file
        """
        # Load files
        msp_df = self.load_data(msp_file)
        # Restrict df to detected genes
        detected_genes = self.gene_count.loc[self.gene_count["value"] > 0, "gene_id"]
        msp_df = msp_df.loc[msp_df["gene_id"].isin(detected_genes)]
        # Restrict df to detected msp
        msp_df = msp_df.loc[
            msp_df["msp_name"].isin(
                self.msp_table.loc[self.msp_table["value"] > 0, "msp_name"]
            )
        ]
        # Merge each provided db
        annot_df = pd.concat(
            [
                self.load_data(db)[["gene_id", "annotation"]]
                for db in annot_file.values()
            ],
            ignore_index=True,
        )
        annot_df = annot_df.loc[annot_df["gene_id"].isin(detected_genes)]
        annotated_msp_df = msp_df.merge(annot_df)
        return annotated_msp_df

    def compute_completeness(
        self, mod: list[set[str]], annotated_msp: pd.DataFrame
    ) -> dict[str, float]:
        """Compute completeness of a given module in all available MSP.

        :param mod: a module (defined by all its alternatives)
        :param annotated_msp: a dataframe containing functionnaly annotated gene content
        of all MSP.
        """
        return (
            annotated_msp.groupby("msp_name")["annotation"]
            .apply(lambda x: self.compute_max(mod, set(x)))
            .to_dict()
        )

    def compute_max(self, mod: list[set[str]], ko: set) -> float:
        """Compute maximum completeness of a module across
        all its alternative according to a set of KO.

        :param mod: a module defined by all its alternatives
        :param ko: a set of ko were the alternatives should be searched for
        """
        return max((len(alt.intersection(ko)) / len(alt) for alt in mod))

    def compute_completeness_all(
        self, all_mod: dict[str, list[set[str]]], annotated_msp: pd.DataFrame
    ) -> dict[str, dict[str, float]]:
        """Compute completeness of all modules in all MSP.

        :param all_mod: dictionnary where key = module id and value = all alternatives
        of the module
        :param annotated_msp: a dataframe with functionnaly annotated genes content of
        all MSP.
        """
        return {
            mod: self.compute_completeness(alt, annotated_msp)
            for (mod, alt) in all_mod.items()
        }

    def compute_module_abundance(
        self,
        msp_file: Path,
        annot_file: dict[str, Path],
        all_mod: dict[str, list[set[str]]],
        completeness: float,
    ) -> None:
        """Compute all modules abundance in the sample.

        :param msp_file: path to the msp definition file
        :param annot_file: path to the gene annotation file (functional)
        :param all_mod: dictionnary of all modules with their list of alternatives
        :param completeness: KO ratio above which a module is set as present in an MSP
        """
        # Merge the data
        annotated_msp = self.merge_catalogue_info(
            msp_file=msp_file, annot_file=annot_file
        )
        # Compute all completeness for all modules and all MSP
        cpltd_dict = self.compute_completeness_all(
            all_mod=all_mod, annotated_msp=annotated_msp
        )
        # Reformat module completeness as dataframe
        self.mod_completeness = (
            pd.DataFrame.from_dict(cpltd_dict)
            .melt(ignore_index=False)
            .reset_index()
            .rename(
                columns={
                    "index": "msp_name",
                    "variable": "mod_id",
                }
            )
        )
        # Restrict to msp whose completeness is above threshold
        mod_dict = {
            mod: {msp for (msp, cmpltd) in msp_dict.items() if cmpltd >= completeness}
            for (mod, msp_dict) in cpltd_dict.items()
        }
        # Compute module abundance
        module_abundance = {
            mod: self.msp_table.loc[
                self.msp_table["msp_name"].isin(msp_set), "value"
            ].sum()
            for (mod, msp_set) in mod_dict.items()
        }
        self.mod_table = (
            pd.DataFrame.from_dict(module_abundance, orient="index", columns=["value"])
            .reset_index()
            .rename(columns={"index": "mod_id"})
        )

    def execute(self) -> None:
        "Normalize the samples and compute MSP and functions abundances."
        # Initialize dictionnary for json file
        config_param = {}
        config_stats = {}
        config_mapping = {"database_type": self.database_type}
        # Part 1: NORMALIZATION
        if self.rarefaction_level > Profiler.NO_RAREFACTION:
            logging.info("Run rarefaction.")
            self.rarefy(
                rarefaction_level=self.rarefaction_level,
                unmapped_reads=self.unmapped_reads,
                seed=self.seed,
            )
        else:
            logging.info("No rarefaction.")
        if self.normalization == "coverage":
            logging.info("Run coverage normalization.")
            self.normalize_coverage(
                trim_length=int(self.sample_config["mapping"]["trim"])
            )
        elif self.normalization == "fpkm":
            logging.info("Run fpkm normalization.")
            self.normalize_fpkm(
                rarefaction_level=self.rarefaction_level,
                unmapped_reads=self.unmapped_reads,
            )
        else:
            logging.info("No normalization.")
        # Write the normalized count table
        logging.info("Save gene table.")
        gene_table_file = self.stage2_dir / f"{self.output_base_filename}_genes.tsv.xz"
        with lzma.open(gene_table_file, "wt", preset=0) as out:
            self.gene_count.to_csv(out, sep="\t", index=False)
        # Update config dictionnary
        config_param.update(
            {
                "normalization": str(self.normalization),
                "rarefaction_level": str(self.rarefaction_level),
                "seed": str(self.seed),
            }
        )
        config_stats["gene_count"] = str(
            len(
                self.gene_count.loc[
                    self.gene_count["value"] > 0,
                    "gene_id",
                ]
            )
        )
        # Part 2: TAXONOMIC PROFILING
        # Restrict to MSP of interest
        logging.info("Get MSP core genes.")
        msp_set = self.get_msp_core(self.msp_filename, self.core_size)
        # Compute MSP
        logging.info("Compute MSP profiles.")
        self.compute_msp(msp_dict=msp_set, filter_pc=self.msp_filter)
        # Write the MSP table
        logging.info("Save MSP profiles.")
        msp_table_file = self.stage2_dir / f"{self.output_base_filename}_msp.tsv.xz"
        with lzma.open(msp_table_file, "wt", preset=0) as out:
            self.msp_table.to_csv(out, sep="\t", index=False)
        # Compute MSP stats
        logging.info("Compute MSP stats.")
        msp_stats = self.compute_msp_stats(self.msp_filename)
        # Update and save config file
        config_param.update(
            {
                "msp_core_size": str(self.core_size),
                "msp_filter": str(self.msp_filter),
                "msp_def": self.msp_filename.name,
            }
        )
        config_stats.update(
            {
                "msp_count": str(
                    len(
                        self.msp_table.loc[
                            self.msp_table["value"] > 0,
                            "msp_name",
                        ]
                    )
                ),
                "msp_signal": str(msp_stats),
            }
        )
        # Part 3: FUNCTIONAL PROFILING
        if self.database_type == "complete":
            single_fun_db = ["mustard", "kegg", "dbcan", "resfinder", "resfinderfg"]
            single_fun_by_msp_db = [
                "mustard",
                "kegg",
                "dbcan",
                "resfinder",
                "resfinderfg",
            ]
            for db, db_filename in self.db_filenames.items():
                # By sum of genes
                if db in single_fun_db:
                    logging.info("Compute %s abundances as sum of gene abundances.", db)
                    self.compute_ko_abundance(annot_file=db_filename)
                    logging.info("Save functional profiles table.")
                    fun_table_file = (
                        self.stage2_dir
                        / f"{self.output_base_filename}_{db}_as_genes_sum.tsv.xz"
                    )
                    with lzma.open(fun_table_file, "wt", preset=0) as out:
                        self.functions.to_csv(out, sep="\t", index=False)
                    # Compute functional statistics
                    functional_stats = self.compute_ko_stats(
                        annot_file=db_filename,
                        by_msp=False,
                        msp_def_filename=self.msp_filename,
                    )

                    # Update config file
                    config_param[f"{db}_filename"] = self.db_filenames[db].name
                    config_stats[f"{db}_signal_by_genes"] = str(functional_stats)
                # By sum of MSPs
                if db in single_fun_by_msp_db:
                    logging.info("Compute %s abundances as sum of MSP abundances.", db)
                    self.compute_ko_abundance_by_msp(
                        annot_file=db_filename, msp_def_filename=self.msp_filename
                    )
                    logging.info("Save functional profiles table.")
                    fun_table_file = (
                        self.stage2_dir
                        / f"{self.output_base_filename}_{db}_as_msp_sum.tsv.xz"
                    )
                    with lzma.open(fun_table_file, "wt", preset=0) as out:
                        self.functions.to_csv(out, sep="\t", index=False)

                    # Compute functinal statistics
                    functional_stats = self.compute_ko_stats(
                        annot_file=db_filename,
                        by_msp=True,
                        msp_def_filename=self.msp_filename,
                    )

                    # Update config file
                    config_param[f"{db}_filename"] = self.db_filenames[db].name
                    config_stats[f"{db}_signal_by_msp"] = str(functional_stats)

            # Part 4 Module computation
            # Get db filenames required for module computation
            module_db_filenames = {
                db: self.db_filenames[db] for db in ["kegg", "eggnog", "tigrfam"]
            }
            logging.info("Parse module definition file.")
            module_dict = Parser(self.module_path)
            module_dict.execute()
            logging.info("Compute modules.")
            self.compute_module_abundance(
                msp_file=self.msp_filename,
                annot_file=module_db_filenames,
                all_mod=module_dict.module_dict_alt,
                completeness=self.completeness,
            )
            module_table_file = (
                self.stage2_dir / f"{self.output_base_filename}_modules.tsv.xz"
            )
            with lzma.open(module_table_file, "wt", preset=0) as out:
                self.mod_table.to_csv(out, sep="\t", index=False)
            module_completeness_file = (
                self.stage2_dir
                / f"{self.output_base_filename}_modules_completeness.tsv.xz"
            )
            with lzma.open(module_completeness_file, "wt", preset=0) as out:
                self.mod_completeness.to_csv(out, sep="\t", index=False)
            # Update config files
            config_param["modules_db"] = ",".join(module_db_filenames.keys())
            config_param["modules_db_filenames"] = ",".join(
                [value.name for value in module_db_filenames.values()]
            )
            config_param["modules_def"] = self.module_path.name
            config_param["module_completeness"] = str(self.completeness)

        # Update and write ini file
        update_config = self.update_json(
            self.sample_config, "profiling_parameters", config_param
        )
        update_config = self.update_json(update_config, "profiling_stats", config_stats)
        update_config = self.update_json(update_config, "mapping", config_mapping)
        self.save_config(
            update_config,
            self.stage2_dir / f"{self.output_base_filename}_census_stage_2.json",
        )

        logging.info("Process ended without errors.")
