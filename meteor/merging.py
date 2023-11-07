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
from configparser import ConfigParser
from dataclasses import dataclass
from typing import Type
import pandas as pd
from pathlib import Path
import logging
import sys
from functools import reduce


@dataclass
class Merging(Session):
    """
    Merge profiled outputs.
    """
    meteor: Type[Component]
    pattern: str
    check_param: bool
    output: Path

    def find_files_to_merge(self, input_dir: Path, pattern: str) -> list[Path]:
        files_to_merge = list(input_dir.glob("*" + pattern))
        return files_to_merge

    def find_associated_ini(self, list_files: list[Path]) -> dict[Path, ConfigParser]:
        return {x: self.read_ini(x.with_suffix(".ini")) for x in list_files}

    def get_sample_name(self, dict_config: dict[Path, ConfigParser]) -> dict[str, Path]:
        return {my_config["sample_info"]["sample_name"]: my_path for my_path, my_config in dict_config.items()}

    def extract_ini_info(self, list_config: dict[Path, ConfigParser], section: str) -> dict[Path, dict[str, str]]:
        "Get all information from the required section"
        try:
            assert all((section in my_config for my_config in list_config.values()))
        except AssertionError:
            logging.error("Missing the section %s in one of the ini file.", section)
            sys.exit()
        return {my_path: dict(my_config[section]) for (my_path, my_config) in list_config.items()}

    def compare_section_info(self, all_sections_info: dict[Path, dict[str, str]]) -> bool:
        "Check that all dictionnary are equal"
        # Define an ini file of reference
        ref_config = list(all_sections_info.keys())[0]
        logging.info("Use the ini file associated to %s as a reference.", ref_config)
        # Compare all other entries to this reference
        ref_comparison_bool = {my_path: my_config == all_sections_info[ref_config]
                               for (my_path, my_config) in all_sections_info.items()}
        # Get ini file for which this is not true
        problematic_config = {key: value for (key, value) in ref_comparison_bool.items() if value is False}
        if len(problematic_config) > 0:
            for i in problematic_config:
                logging.error("The file %s does not match the reference ini file.", i)
            sys.exit()
        else:
            logging.info("All files have the same parameters.")
            return True

    def merge_df(self, dict_path: dict[str, Path], key_merging: str) -> pd.DataFrame:
        "Merge many data frames contained in a list"
        # Load the data frames
        list_df = [pd.read_table(my_path).rename(columns={
            self.meteor.value_column: my_sample})[[key_merging, my_sample]]
            for (my_sample, my_path) in dict_path.items()]
        merged_df = reduce(lambda left, right: pd.merge(left,right,on=key_merging,how="outer"), list_df)
        return merged_df

    def execute(self):
        # Fetch the data
        logging.info("Fetching files matching the pattern %s...", self.pattern)
        files_to_merge = self.find_files_to_merge(input_dir = self.meteor.mapping_dir, pattern = self.pattern)
        logging.info("There was %s files that correspond to the pattern.", len(files_to_merge))
        # Get the associated ini file (obligatory to get the sample names)
        all_config = self.find_associated_ini(files_to_merge)
        # Check the parameters if required
        if self.check_param:
            # Check the profiling parameters section
            extracted_info = self.extract_ini_info(all_config, "profiling_parameters")
            self.compare_section_info(extracted_info)
        else:
            logging.info("Their parameters will not be checked prior to merging.")
        # Guess the key that should be used
        logging.info("Look for the common key for merging.")
        with open(files_to_merge[0], encoding="UTF-8") as ref:
            real_colnames = set(ref.readline().strip("\n").split("\t"))
        possible_keys = [self.meteor.msp_column, self.meteor.gene_column,
                         self.meteor.ko_column, self.meteor.module_column]
        key_merging_set = set(real_colnames).intersection(possible_keys)
        try:
            assert len(key_merging_set) == 1
        except AssertionError:
            logging.error("No unique key was found in the columns you try to merge: %s", ", ".join(real_colnames))
            sys.exit()
        key_merging = list(key_merging_set)[0]
        logging.info("Key that will be used for merging: %s", key_merging)
        # Get sample names
        logging.info("Fetching sample names.")
        sample_names = self.get_sample_name(all_config)
        # Perform the actual merging
        logging.info("Merging the data.")
        merged_df = self.merge_df(sample_names, key_merging)
        # Save the data frame
        logging.info("Save the data.")
        merged_df.to_csv(self.output, sep = "\t", index = False)




