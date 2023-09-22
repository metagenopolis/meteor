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
from dataclasses import dataclass
from typing import Type
import pandas as pd
from pathlib import Path
import numpy as np
import logging
import os
import sys
from datetime import date


@dataclass
class Profiler(Session):
    """Profile session for abundance and annotation
    """
    meteor: Type[Component]
    input_count_table: Path
    input_ini: Path
    suffix_file: str
    rarefaction_level: int
    seed: int
    normalization: str
    compute_msp_bool: bool
    core_size: int
    msp_filter: float
    compute_functions_bool: bool
    annot_db: str
    by_msp: bool
    compute_modules_bool: bool
    module_path: Path
    module_db: str
    completeness: float

    def __post_init__(self):
        # Check the input count table
        self.check_file(self.input_count_table, {self.meteor.gene_column,
                                                 self.meteor.value_column,
                                                 self.meteor.gene_length_column})
        # Load the count table
        self.gene_count = pd.read_table(self.input_count_table)
        self.count_column = self.gene_count.columns[2]
        self.gene_count[self.meteor.value_column] = self.gene_count[self.meteor.value_column].round(0).astype("int")

        # Initialize the ini file
        if self.input_ini is None:
            self.input_ini = Path(self.input_count_table).with_suffix(".ini")
        self.sample_config = self.read_ini(self.input_ini)

        # Add session info
        config_session = {}
        config_session["user"] = os.getlogin()
        config_session["date"] = str(date.today())
        self.sample_config = self.update_ini(self.sample_config, "profiling_session", config_session)

        # Initialize the unmapped_count
        try:
            self.unmapped_reads = int(self.sample_config["mapping"]["unmapped_reads"])
        except KeyError:
            logging.info("No unmapped reads information found. Set to 0.")
            self.unmapped_reads = 0

        # Initialize the sample name
        self.sample_name = self.sample_config["sample_info"]["sample_name"]

        # Initialize the ini ref config parser:
        self.ref_config = self.get_reference_info(self.meteor.ref_dir)

        # Define MSP filename
        if self.by_msp or self.compute_msp_bool:
            try:
                self.msp_filename = (self.meteor.ref_dir /
                                     self.ref_config["reference_file"]["database_dir"] /
                                     self.ref_config[self.meteor.ko_column]["msp"])
                self.msp_filename = Path(self.msp_filename)
                assert self.msp_filename.is_file()
            except KeyError:
                logging.error("No MSP file provided in the reference ini file.")
                sys.exit()
            except AssertionError:
                logging.error("The MSP file %s does not exist.", self.msp_filename)
                sys.exit()
            self.check_file(self.msp_filename, {self.meteor.msp_column,
                                                self.meteor.gene_column,
                                                self.meteor.gene_class_column})
        else:
            self.msp_filename = None

        # List the functional db
        if self.compute_functions_bool:
            self.annot_db_list = self.annot_db.split(",")
        else:
            self.annot_db_list = []
        if self.compute_modules_bool:
            self.module_db_list = self.module_db.split(",")
        else:
            self.module_db_list = []
        try:
            self.annot_db_dict = {db: Path(self.meteor.ref_dir /
                                        self.ref_config["reference_file"]["database_dir"] /
                                        self.ref_config[self.meteor.ko_column][db])
                                for db
                                in self.annot_db_list}
            self.module_db_dict = {db: Path(self.meteor.ref_dir /
                                            self.ref_config["reference_file"]["database_dir"] /
                                            self.ref_config[self.meteor.ko_column][db])
                                for db
                                in self.module_db_list}
            assert all(x.is_file() for x in self.annot_db_dict.values())
            assert all(x.is_file() for x in self.module_db_dict.values())
        except KeyError:
            logging.error("Missing annotation databases in the reference ini file.")
            sys.exit()
        except AssertionError:
            logging.error("The annotation files does not exist.")
            sys.exit()
        for db in list(self.annot_db_dict.values()) + list(self.module_db_dict.values()):
            self.check_file(db, {self.meteor.gene_column, self.meteor.ko_column})

        # Define output names
        gene_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_norm.tsv"
        msp_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_msp.tsv"
        functions_table_output = {db: self.meteor.mapping_dir /
                                  f"{self.sample_name}_{self.suffix_file}_{db}_functions.tsv"
                                  for db
                                  in self.annot_db_list}
        modules_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_modules.tsv"
        self.output_filenames = {}
        self.output_filenames["gene_table_norm"] = gene_table_output
        self.output_filenames["msp_table"] = msp_table_output
        self.output_filenames["functions_table"] = functions_table_output
        self.output_filenames["modules_table"] = modules_table_output

        # Initialize the module definition file
        if self.compute_modules_bool:
            if self.module_path is None:
                self.module_path = Path(os.path.join(os.path.dirname(__file__),
                                                    "all_modules_definition_GMM_GBM_KEGG_107.tsv"))
            try:
                assert self.module_path.is_file()
            except AssertionError:
                logging.error("The file %s does not exist.", self.module_path)
            # Do not check columns here because no header

    def rarefy(self, rarefaction_level: int, unmapped_reads: int, seed: int) -> None:
        count_column = self.meteor.value_column
        # Add the unmapped count
        self.gene_count.loc[len(self.gene_count)] = {self.meteor.gene_column: -1,
                                                     self.meteor.gene_length_column: 1000,
                                                     count_column: unmapped_reads}

        # Check if rarefaction is possible
        if self.gene_count[count_column].sum() > rarefaction_level:
            # Transform into a long array of gene_id according to gene occurence
            array_to_rarefy = np.repeat(self.gene_count[self.meteor.gene_column],
                                        self.gene_count[count_column])
            # Randomly choose among the long array the selected reads
            rng = np.random.default_rng(seed = seed)
            array_rarefied = rng.choice(array_to_rarefy,
                                        size = rarefaction_level,
                                        replace = False)
            # Count gene_id occurence to get back to short gene_id list
            unique, counts = np.unique(array_rarefied, return_counts=True)
            self.gene_count[count_column] = np.zeros(self.gene_count.shape[0], dtype = "int")
            self.gene_count.loc[np.in1d(self.gene_count[self.meteor.gene_column], unique), count_column] = counts

        # Remove the counts for the gene "-1" (unmapped_reads)
        self.gene_count = self.gene_count[self.gene_count[self.meteor.gene_column] != -1]

    def normalize_coverage(self) -> None:
        """Normalize by coverage
        """
        count_column = self.meteor.value_column
        self.gene_count[count_column] = (self.gene_count[count_column] /
                                         self.gene_count[self.meteor.gene_length_column] * 1000.)

    def normalize_fpkm(self, rarefaction_level: int, unmapped_reads: int) -> None:
        """Normalize matrix using fpkm method

        :param rarefaction_level: Value of rarefaction level
        :param unmapped_reads: Value of the unmapped reads
        """
        count_column = self.meteor.value_column
        # Compute the unmapped reads after rarefaction
        if rarefaction_level > 0:
            unmapped_reads_after_rf = self.rarefaction_level - self.gene_count[count_column].sum()
        else:
            unmapped_reads_after_rf = unmapped_reads
        ### Add the unmapped reads after rf to the gene count table as a pseudo gene
        self.gene_count.loc[len(self.gene_count)] = {self.meteor.gene_column: -1,
                                                     self.meteor.gene_length_column: 1000,
                                                     count_column: unmapped_reads_after_rf}
        ### Normalize
        self.gene_count[count_column] = self.gene_count[count_column] / self.gene_count[self.meteor.gene_length_column]
        self.gene_count[count_column] = self.gene_count[count_column].div(self.gene_count[count_column].sum())

    def compute_msp(self, msp_dict: dict[str, set[str]], filter_pc: float) -> None:
        count_column = self.meteor.value_column
        # Compute how many genes are seen for a given msp
        # Restrict to gene table to core genes
        all_core_genes = {item for sublist in msp_dict.values() for item in sublist}
        gene_count_core = self.gene_count.loc[self.gene_count[self.meteor.gene_column].isin(all_core_genes)]
        msp_filter = {msp:(gene_count_core.loc[gene_count_core[self.meteor.gene_column].isin(set_genes), count_column] > 0).sum() /
                      len(set_genes)
                      for (msp, set_genes)
                      in msp_dict.items()}
        # Compute mean abundance if gene count is above filter threshold, otherwise 0
        msp_table_dict = {
            msp:(gene_count_core.loc[gene_count_core[self.meteor.gene_column].isin(set_genes), count_column].mean()
                 if msp_filter[msp] >= filter_pc else 0)
                 for (msp, set_genes) in msp_dict.items()}
        self.msp_table = pd.DataFrame.from_dict(msp_table_dict,
                                                orient = "index",
                                                columns = [self.meteor.value_column]).reset_index().rename(columns={"index": self.meteor.msp_column})

    def get_msp_core(self, msp_def_filename: Path, core_size: int) -> dict:
        # Load msp file
        msp_df = pd.read_table(msp_def_filename)
        # Restrict to core
        msp_df_selection = msp_df.loc[msp_df[self.meteor.gene_class_column] == "core"]
        # Return the df as a dict of set
        msp_dict = msp_df_selection.groupby(self.meteor.msp_column)[self.meteor.gene_column].apply(lambda x: set(x.head(core_size))).to_dict()
        return msp_dict

    def compute_msp_stats(self, msp_def_filename: Path) -> float:
        count_column = self.meteor.value_column
        # Load msp file
        msp_df = pd.read_table(msp_def_filename)
        # Get the ensemble of genes used in MSP
        all_msp_genes = msp_df[self.meteor.gene_column].unique()
        # Get the percentage of reads that map on an MSP
        msp_reads_pc = (self.gene_count.loc[self.gene_count[self.meteor.gene_column].isin(all_msp_genes), count_column].sum() /
                        self.gene_count[count_column].sum())
        return round(msp_reads_pc, 2)

    def compute_ko_abundance(self, annot_file: Path) -> None:
        count_column = self.meteor.value_column
        # Load annotation file
        annot_df = pd.read_table(annot_file)
        # Merge count table and gene annotation
        merged_df = pd.merge(annot_df, self.gene_count, left_on=self.meteor.gene_column, right_on=self.meteor.gene_column)
        # Compute sum of KO
        aggregated_count = merged_df.groupby(self.meteor.ko_column)[count_column].sum().reset_index()
        self.functions = aggregated_count

    def compute_ko_abundance_by_msp(self, annot_file: Path, msp_def_filename: Path) -> None:
        # Load annotation file
        annot_df = pd.read_table(annot_file)
        # Load MSP file
        msp_df = pd.read_table(msp_def_filename)
        # Merge both data frames
        msp_df_annotated = pd.merge(msp_df, annot_df)
        # Create a dict ko: {msp1, msp2}
        ko_dict = msp_df_annotated.groupby(self.meteor.ko_column)[self.meteor.msp_column].apply(set).to_dict()
        # Compute abundance based on MSP abundance
        ko_dict_ab = {ko: self.msp_table.loc[self.msp_table[self.meteor.msp_column].isin(msp_set), self.meteor.value_column].sum()
                      for (ko, msp_set)
                      in ko_dict.items()}
        self.functions = pd.DataFrame.from_dict(ko_dict_ab, orient = "index", columns = [self.meteor.value_column]).reset_index().rename(columns={"index": self.meteor.ko_column})

    def compute_ko_stats(self, annot_file: Path, by_msp: bool, msp_def_filename: Path) -> float:
        count_column = self.meteor.value_column
        # Load annotation file
        annot_df = pd.read_table(annot_file)
        if by_msp:
            # Load MSP file
            msp_df = pd.read_table(msp_def_filename)
            # Merge both data frames
            annot_df = pd.merge(msp_df, annot_df)
        # Get the genes in MSP AND annotated
        all_msp_genes = annot_df[self.meteor.gene_column].unique()
        # Get the percentage of reads that map on these genes
        annot_reads_pc = (self.gene_count.loc[self.gene_count[self.meteor.gene_column].isin(all_msp_genes), count_column].sum() /
                          self.gene_count[count_column].sum())
        return round(annot_reads_pc, 2)

    def merge_catalogue_info(self, msp_file: Path, annot_file: dict[str, Path]) -> pd.DataFrame:
        count_column = self.meteor.value_column
        # Load files
        msp_df = pd.read_table(msp_file)
        # Restrict df to detected genes
        detected_genes = self.gene_count.loc[self.gene_count[count_column] > 0, self.meteor.gene_column]
        msp_df = msp_df.loc[msp_df[self.meteor.gene_column].isin(detected_genes)]
        # Restrict df to detected msp
        msp_df = msp_df.loc[msp_df[self.meteor.msp_column].isin(self.msp_table.loc[self.msp_table[self.meteor.value_column] > 0, self.meteor.msp_column])]
        # Merge each provided db
        annot_df = pd.concat([pd.read_table(db)[[self.meteor.gene_column, self.meteor.ko_column]]
                              for db
                              in annot_file.values()],
                             ignore_index = True)
        annot_df = annot_df.loc[annot_df[self.meteor.gene_column].isin(detected_genes)]
        annotated_msp_df = msp_df.merge(annot_df)
        return annotated_msp_df

    def compute_completeness(self, mod: list[set[str]], annotated_msp: pd.DataFrame) -> dict[str, float]:
        "Compute completeness of a given module in all available MSP."
        return annotated_msp.groupby(self.meteor.msp_column)[self.meteor.ko_column].apply(lambda x: self.compute_max(mod, set(x))).to_dict()

    def compute_max(self, mod: list[set[str]], ko: set) -> float:
        "Compute maximum completeness of a module across all its alternative according to a set of KO."
        return max((len(alt.intersection(ko)) / len(alt) for alt in mod))

    def compute_completeness_all(self,
                                 all_mod: dict[str, list[set[str]]],
                                 annotated_msp: pd.DataFrame) -> dict[str, dict[str, float]]:
        "Compute completeness of all modules in all MSP."
        return {mod: self.compute_completeness(alt, annotated_msp) for (mod, alt) in all_mod.items()}

    def compute_module_abundance(self, msp_file: Path, annot_file: dict[str, Path],
                                 all_mod: dict[str, list[set[str]]],
                                 completeness: float) -> None:
        "Compute all modules abundance in the sample."
        count_column = self.meteor.value_column
        # Merge the data
        annotated_msp = self.merge_catalogue_info(msp_file=msp_file, annot_file=annot_file)
        # Compute all completeness for all modules and all MSP
        cpltd_dict = self.compute_completeness_all(all_mod=all_mod, annotated_msp=annotated_msp)
        # Restrict to msp whose completeness is above threshold
        mod_dict = {mod: {msp for (msp, cmpltd) in msp_dict.items() if cmpltd >= completeness}
                    for (mod, msp_dict)
                    in cpltd_dict.items()}
        # Compute module abundance
        module_abundance = {mod: self.msp_table.loc[self.msp_table[self.meteor.msp_column].isin(msp_set), self.meteor.value_column].sum()
                            for (mod, msp_set)
                            in mod_dict.items()}
        self.mod_table = pd.DataFrame.from_dict(module_abundance, orient = "index", columns = [self.meteor.value_column]).reset_index().rename(columns={"index": self.meteor.module_column})


    def execute(self) -> bool:
        "Normalize the samples and compute MSP and functions abundances."
        # Part 1: NORMALIZATION
        if self.rarefaction_level > 0:
            logging.info("Run rarefaction.")
            self.rarefy(rarefaction_level=self.rarefaction_level,
                        unmapped_reads=self.unmapped_reads,
                        seed=self.seed)
        else:
            logging.info("No rarefaction.")
        if self.normalization == "coverage":
            logging.info("Run coverage normalization.")
            self.normalize_coverage()
        elif self.normalization == "fpkm":
            logging.info("Run fpkm normalization.")
            self.normalize_fpkm(rarefaction_level=self.rarefaction_level,
                                unmapped_reads=self.unmapped_reads)
        else:
            logging.info("No normalization.")
        # Write the normalized count table
        logging.info("Save gene table.")
        self.gene_count.to_csv(self.output_filenames["gene_table_norm"], sep = "\t", index = False)
        # Update config file
        config_norm = {}
        config_norm["normalization"] = self.normalization
        config_norm["rarefaction_level"] = str(self.rarefaction_level)
        config_norm["seed"] = str(self.seed)
        self.sample_config = self.update_ini(self.sample_config, "profiling_parameters", config_norm)
        self.save_config(self.sample_config, Path(self.output_filenames["gene_table_norm"]).with_suffix(".ini"))

        # Part 2: TAXONOMIC PROFILING
        if self.compute_msp_bool or self.by_msp:
            if not self.compute_msp_bool:
                logging.info("MSP abundances will be computed since you requested to compute functions via MSP.")
            # Restrict to MSP of interest
            logging.info("Get MSP core genes.")
            msp_set = self.get_msp_core(self.msp_filename, self.core_size)
            # Compute MSP
            logging.info("Compute MSP profiles.")
            self.compute_msp(msp_dict=msp_set, filter_pc=self.msp_filter)
            # Write the MSP table
            logging.info("Save MSP profiles.")
            self.msp_table.to_csv(self.output_filenames["msp_table"], sep = "\t", index = False)
            # Compute MSP stats
            logging.info("Compute MSP stats.")
            msp_stats = self.compute_msp_stats(self.msp_filename)
            # Update and save config file
            config_param_msp = {}
            config_param_msp["msp_core_size"] = str(self.core_size)
            config_param_msp["msp_filter"] = str(self.msp_filter)
            config_param_msp["msp_def"] = self.msp_filename.name
            config_stats_msp = {}
            config_stats_msp["msp_count"] = str(len(self.msp_table.loc[self.msp_table[self.meteor.value_column] > 0, self.meteor.msp_column]))
            config_stats_msp["msp_signal"] = str(msp_stats)
            self.sample_config = self.update_ini(self.sample_config, "profiling_parameters", config_param_msp)
            self.sample_config = self.update_ini(self.sample_config, "profiling_stats", config_stats_msp)
            self.save_config(self.sample_config, Path(self.output_filenames["msp_table"]).with_suffix(".ini"))

        # Part 3: FUNCTIONAL PROFILING
        if self.compute_functions_bool:
            # Get functions to use
            for (db, annot_file) in self.annot_db_dict.items():
                logging.info("Compute %s abundances.", db)
                if self.by_msp:
                    logging.info("Use MSP abundances. Remove --by_msp not to use them.")
                    self.compute_ko_abundance_by_msp(annot_file=annot_file, msp_def_filename=self.msp_filename)
                else:
                    logging.info("MSP abundance are not used. Add --by_msp to use them for computations.")
                    self.compute_ko_abundance(annot_file=annot_file)
                logging.info("Save annotation file.")
                self.functions.to_csv(self.output_filenames["functions_table"][db], sep = "\t", index=False)
                # Compute functinal statistics
                functional_stats = self.compute_ko_stats(annot_file=annot_file,
                                                         by_msp=self.by_msp, msp_def_filename=self.msp_filename)
                # Update config file
                config_param_db = {}
                config_param_db["function_db"] = db
                config_param_db["function_filename"] = annot_file.name
                config_param_db["function_by_msp"] = str(self.by_msp)
                config_stats_db = {}
                config_stats_db["function_signal"] = str(functional_stats)
                update_config = self.update_ini(self.sample_config, "profiling_parameters", config_param_db)
                update_config = self.update_ini(update_config, "profiling_stats", config_stats_db)
                self.save_config(update_config,
                                 Path(self.output_filenames["functions_table"][db]).with_suffix(".ini"))

        # Part 4 Module computation
        if self.compute_modules_bool:
            logging.info("Compute module abundances using the file %s", self.module_path.name)
            logging.info("Parse module definition file.")
            module_dict = Parser(self.module_path)
            module_dict.execute()
            logging.info("Compute modules.")
            self.compute_module_abundance(msp_file=self.msp_filename,
                                          annot_file=self.module_db_dict,
                                          all_mod=module_dict.module_dict_alt,
                                          completeness=self.completeness)
            self.mod_table.to_csv(self.output_filenames["modules_table"], sep = "\t", index = False)
            # Update config files
            config_param_modules = {}
            config_param_modules["modules_def"] = self.module_path.name
            config_param_modules["module_completeness"] = str(self.completeness)
            update_config = self.update_ini(self.sample_config, "profiling_parameters", config_param_modules)
            self.save_config(update_config,
                                 Path(self.output_filenames["modules_table"]).with_suffix(".ini"))

        logging.info("Process ended without errors.")

        return True
