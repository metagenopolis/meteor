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
from pathlib import Path
import numpy as np
import logging
from os import sys, getlogin
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
    normalization: str
    compute_mgs_bool: bool
    core_size: int
    mgs_filter: float
    compute_functions_bool: bool
    annot_db: str
    by_mgs: bool
    compute_modules_bool: bool
    module_path: Path
    module_db: str
    completude: float

    def __post_init__(self):
        # Load the count table
        self.gene_count = pd.read_table(self.input_count_table)
        self.count_column = self.gene_count.columns[2]
        self.gene_count[self.count_column] = self.gene_count[self.count_column].round(0).astype("int")

        # Initialize the ini file
        if self.input_ini is None:
            self.input_ini = Path(self.input_count_table).with_suffix(".ini")
        self.sample_config = self.read_ini(self.input_ini)
            
        # Initialize the unmapped_count
        try:
            self.unmapped_reads = int(self.sample_config["mapping"]["unmapped_reads"])
        except KeyError:
            logging.info("No unmapped reads information found. Set to 0.")
            self.unmapped_reads = 0

        # Initialize the sample name
        self.sample_name = self.sample_config["sample_info"]["sample_name"]

        # List the annot_db
        self.annot_db_list = self.annot_db.split(",")

        # Define output names
        gene_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_norm.tsv"
        mgs_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_mgs.tsv"
        functions_table_output = {db: self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_{db}_functions.tsv" for db in self.annot_db_list}
        modules_table_output = self.meteor.mapping_dir / f"{self.sample_name}_{self.suffix_file}_modules.tsv"
        self.output_filenames = {}
        self.output_filenames["gene_table_norm"] = gene_table_output
        self.output_filenames["mgs_table"] = mgs_table_output
        self.output_filenames["functions_table"] = functions_table_output
        self.output_filenames["modules_table"] = modules_table_output

        # Initialize the ini ref config parser:
        self.ref_config = self.get_reference_info(self.meteor.ref_dir)


    def rarefy(self, count_column: str, rarefaction_level: int, unmapped_reads: int) -> None:
        # Add the unmapped count
        self.gene_count.loc[len(self.gene_count)] = {'genes_id': -1, 
                                                     'gene_size': 1000, 
                                                     count_column: unmapped_reads}
        
        # Check if rarefaction is possible
        if (self.gene_count[count_column].sum() > rarefaction_level):
            # Transform into a long array of gene_id according to gene occurence
            array_to_rarefy = np.repeat(self.gene_count["genes_id"],
                                        self.gene_count[count_column])
            # Randomly choose among the long array the selected reads
            rng = np.random.default_rng()
            array_rarefied = rng.choice(array_to_rarefy, 
                                        size = rarefaction_level, 
                                        replace = False)
            # Count gene_id occurence to get back to short gene_id list
            unique, counts = np.unique(array_rarefied, return_counts=True)
            self.gene_count[count_column] = np.zeros(self.gene_count.shape[0], dtype = 'int')
            self.gene_count.loc[np.in1d(self.gene_count['genes_id'], unique), count_column] = counts
            
        # Remove the counts for the gene "-1" (unmapped_reads)
        self.gene_count = self.gene_count[self.gene_count.genes_id != -1]
        
    def normalize_coverage(self, count_column: str) -> None:
        self.gene_count[count_column] = self.gene_count[count_column] / self.gene_count["gene_size"] * 1000

    def normalize_fpkm(self, count_column: str, rarefaction_level: int, unmapped_reads: int) -> None:
        # Compute the unmapped reads after rarefaction
        if rarefaction_level > 0:
            unmapped_reads_after_rf = self.rarefaction_level - self.gene_count[count_column].sum()
        else:
            unmapped_reads_after_rf = unmapped_reads
        ### Add the unmapped reads after rf to the gene count table as a pseudo gene
        self.gene_count.loc[len(self.gene_count)] = {'genes_id': -1, 
                                                     'gene_size': 1000, 
                                                     count_column: unmapped_reads_after_rf}
        ### Normalize
        self.gene_count[count_column] = self.gene_count[count_column] / self.gene_count["gene_size"]
        self.gene_count[count_column] = self.gene_count[count_column].div(self.gene_count[count_column].sum())




    def compute_mgs(self, count_column: str, mgs_dict: dict[str, set[str]], filter: int) -> None:
        # Compute how many genes are seen for a given mgs
        # Restrict to gene table to core genes
        all_core_genes = set([item for sublist in mgs_dict.values() for item in sublist])
        gene_count_core = self.gene_count.loc[self.gene_count["genes_id"].isin(all_core_genes)]
        mgs_filter = {mgs:(gene_count_core.loc[gene_count_core["genes_id"].isin(set_genes), count_column] > 0).sum() / len(set_genes) for (mgs, set_genes) in mgs_dict.items()}
        # Compute mean abundance if gene count is above filter threshold, otherwise 0
        mgs_table_dict = {mgs:(gene_count_core.loc[gene_count_core["genes_id"].isin(set_genes), count_column].mean() if mgs_filter[mgs] >= filter else 0) for (mgs, set_genes) in mgs_dict.items()}
        self.mgs_table = pd.DataFrame.from_dict(mgs_table_dict, orient = "index", columns = ["value"]).reset_index()

    def get_mgs_core(self, mgs_def_filename: Path, core_size: int) -> dict:
        # Load mgs file
        mgs_df = pd.read_table(mgs_def_filename)
        # Restrict to core
        mgs_df_selection = mgs_df.loc[mgs_df["gene_category"] == "core",]
        # Restrict to more connected genes
        mgs_df_selection = mgs_df_selection.groupby("msp_name").head(core_size)
        # Return the df as a dict of set
        return mgs_df_selection.groupby(["msp_name"])["gene_id"].apply(lambda grp: set(grp.value_counts().index)).to_dict()
    
    def compute_mgs_stats(self, count_column: str, mgs_def_filename: Path) -> float:
        # Load mgs file
        mgs_df = pd.read_table(mgs_def_filename)
        # Get the ensemble of genes used in MGS
        all_mgs_genes = mgs_df["gene_id"].unique()
        # Get the percentage of reads that map on an MGS
        mgs_reads_pc = self.gene_count.loc[self.gene_count["genes_id"].isin(all_mgs_genes), count_column].sum() / self.gene_count[count_column].sum()
        return round(mgs_reads_pc, 2)
    
    def compute_ko_abundance(self, count_column: str, annot_file: Path) -> None:
        # Load annotation file
        annot_df = pd.read_table(annot_file)
        # Merge count table and gene annotation
        merged_df = pd.merge(annot_df, self.gene_count, left_on='gene_id', right_on="genes_id")
        # Compute sum of KO
        aggregated_count = merged_df.groupby('annotation')[count_column].sum()
        self.functions = aggregated_count

    

    def execute(self) -> bool:
        "Normalize the samples and compute MGS and functions abundances."
        # Normalize the samples
        if self.rarefaction_level > 0:
            logging.info("Run rarefaction.")
            self.rarefy(count_column=self.count_column,
                        rarefaction_level=self.rarefaction_level,
                        unmapped_reads=self.unmapped_reads)
        else:
            logging.info("No rarefaction.")
        if self.normalization == "coverage":
            logging.info("Run coverage normalization.")
            self.normalize_coverage(count_column=self.count_column)
        elif self.normalization == "fpkm":
            logging.info("Run fpkm normalization.")
            self.normalize_fpkm(count_column=self.count_column, 
                                rarefaction_level=self.rarefaction_level, 
                                unmapped_reads=self.unmapped_reads)
        else:
            logging.info("No normalization.")
        # Write the normalized count table
        logging.info("Save gene table.")
        self.gene_count.to_csv(self.output_filenames["gene_table_norm"], sep = "\t", index = False)
        # Update config file
        config_norm = {}
        config_norm["user"] = getlogin()
        config_norm["date"] = date.today()
        config_norm["normalization"] = self.normalization
        config_norm["rarefaction_level"] = self.rarefaction_level
        self.sample_config = self.update_ini(self.sample_config, "normalization", config_norm)
        self.save_config(self.sample_config, Path(self.output_filenames["gene_table_norm"]).with_suffix(".ini"))
        
        # Compute MGS
        if self.compute_mgs_bool:
            # Get MGS filename
            mgs_filename = self.meteor.ref_dir / self.ref_config["reference_file"]["database_dir"] / self.ref_config["annotation"]["msp"]
            # Restrict to MGS of interest
            logging.info("Get MGS core genes.")
            mgs_set = self.get_mgs_core(mgs_filename, self.core_size)
            # Compute MGS
            logging.info("Compute MGS profiles.")
            self.compute_mgs(count_column=self.count_column, mgs_dict=mgs_set, filter=self.mgs_filter)
            # Write the MGS table
            logging.info("Save MGS profiles.")
            self.mgs_table.to_csv(self.output_filenames["mgs_table"], sep = "\t", index = False)
            # Compute MGS stats
            logging.info("Compute MGS stats.")
            mgs_stats = self.compute_mgs_stats(self.count_column, mgs_filename)
            # Update and save config file
            config_mgs = {}
            config_mgs["user"] = getlogin()
            config_mgs["date"] = date.today()
            config_mgs["mgs_signal"] = mgs_stats
            config_mgs["core_size"] = self.core_size
            config_mgs["filter"] = self.mgs_filter
            self.sample_config = self.update_ini(self.sample_config, "mgs", config_mgs)
            self.save_config(self.sample_config, Path(self.output_filenames["mgs_table"]).with_suffix(".ini"))

        if self.compute_functions_bool:
            # Get functions to use
            for db in self.annot_db_list:
                logging.info(f"Compute {db} abundances.")
                annot_file = self.meteor.ref_dir / self.ref_config["reference_file"]["database_dir"] / self.ref_config["annotation"][db]
                self.compute_ko_abundance(count_column=self.count_column, annot_file=annot_file)
                logging.info("Save annotation file.")
                self.functions.to_csv(self.output_filenames["functions_table"][db], sep = "\t")
                # Update config file
                config_db = {}
                config_db["user"] = getlogin()
                config_db["date"] = date.today()
                config_db["db"] = db
                config_db["filename"] = annot_file
                config_db["by_mgs"] = "no"
                self.sample_config = self.update_ini(self.sample_config, "annotation", config_db)
                self.save_config(self.sample_config, Path(self.output_filenames["functions_table"][db]).with_suffix(".ini"))



        logging.info("Process ended without errors.")
        
        return True
