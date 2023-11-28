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

"""Merge profiles"""

from meteor.session import Session, Component
from configparser import ConfigParser
from dataclasses import dataclass
from typing import Type
import pandas as pd
from pathlib import Path
import logging
import sys
import re
from functools import reduce


@dataclass
class Merging(Session):
    """
    Merge profiled outputs.
    """

    meteor: Type[Component]
    output: Path
    prefix: str
    fast: bool = False

    def __post_init__(self) -> None:
        self.output.mkdir(exist_ok=True)

    def extract_census_stage(self, census_list: list[Path]) -> list[int]:
        """For a list of census files formatted as
        sample1_census_stage_1.ini or sample2_census_stage_2.ini,
        extract the stage number (here [1, 2]).

        :param census_list: list of census_stage.ini files as Path.
        """
        motif = re.compile(r"census_stage_(\d+)\.ini")
        list_match = [motif.search(my_file.name) for my_file in census_list]
        return [int(my_match.group(1)) if my_match else None for my_match in list_match]

    def find_files_to_merge(
        self, input_dir: dict[str, Path], pattern: str
    ) -> dict[str, Path]:
        """Find all files in a list of directories matching a specific pattern.

        :param input_dir: dictionnary of path where files should be searched for
        :param pattern: only filenames matching this pattern will be returned
        """
        dict_to_merge = {
            my_sample: list(my_dir.glob("*" + pattern))
            for my_sample, my_dir in input_dir.items()
        }
        # Check that there is exactly one element in each list
        len_list = list(set([len(value) for value in list(dict_to_merge.values())]))
        assert len(len_list) == 1
        assert len_list[0] == 1
        files_to_merge = {
            my_sample: my_list[0] for my_sample, my_list in dict_to_merge.items()
        }
        return files_to_merge

    def extract_ini_info(
        self, config: ConfigParser, param_dict: dict[str, list[str]]
    ) -> dict[str, str]:
        """Get information from the required section and fields

        :param config: A config parser
        :param param_dict: A dictionnary in the form of {section1: [fields1, fields2], section2: []}
        giving the information that should be extracted from the config.
        """
        # Check that sections are present
        try:
            assert all([my_section in config for my_section in list(param_dict.keys())])
        except AssertionError:
            logging.error("Missing required section in census ini file.")
            sys.exit()
        # If fields is empty, consider all fields of the section
        param_dict = {
            key: list(config[key]) if len(value) == 0 else value
            for key, value in param_dict.items()
        }
        # Check that required fields are present
        try:
            assert all(
                [
                    config.has_option(my_section, my_field)
                    for my_section in param_dict
                    for my_field in param_dict[my_section]
                ]
            )
        except AssertionError:
            logging.error("Missing required fields in census ini file.")
            sys.exit()
        # Get values of each field (remove the section information)
        information = {
            my_field: config[my_section][my_field]
            for my_section in param_dict
            for my_field in param_dict[my_section]
        }
        return information

    def compare_section_info(
        self, all_sections_info: dict[Path, dict[str, str]]
    ) -> int:
        """Check that all provided dictionnary are the same.
        Return the number of config parser that are different
        from the reference one (chosen randomly).

        :param all_sections_info: dictionnary of all config parser information,
        where config parser information is stored as a dictionnary field: value
        """
        # Define an ini file of reference
        ref_config = list(all_sections_info.keys())[0]
        logging.info("Use the ini file associated to %s as a reference.", ref_config)
        # Compare all other entries to this reference
        ref_comparison_bool = {
            my_path: my_config == all_sections_info[ref_config]
            for (my_path, my_config) in all_sections_info.items()
        }
        # Get ini file for which this is not true
        problematic_config = {
            key: value for (key, value) in ref_comparison_bool.items() if not value
        }
        if len(problematic_config) > 0:
            logging.warning(
                "There was %s files that do not match the reference ini file",
                str(len(problematic_config)),
            )
            for i in problematic_config:
                logging.warning("The file %s does not match the reference ini file.", i)
        else:
            logging.info("All files have the same parameters.")
        return len(problematic_config)

    def merge_df(
        self, dict_path: dict[str, Path], key_merging: list[str]
    ) -> pd.DataFrame:
        "Merge many data frames contained in a list"
        """Merge data frames contained in a list according to a pivot key.

        :param dict_path: dictionnary of path of files to merge, where key is sample name
        and will be used as colnames.
        :param key_merging: list of keys that should be used as pivot for merging.
        """
        # Load the data frames
        list_df = [
            pd.read_table(my_path).rename(
                columns={self.meteor.value_column: my_sample}
            )[key_merging + [my_sample]]
            for (my_sample, my_path) in dict_path.items()
        ]
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, how="outer"),
            list_df,
        )
        return merged_df

    def execute(self):
        "Merge all files generated by either mapper or profiler."
        # Fetch all census ini files
        all_census = list(Path(self.meteor.profile_dir).glob("**/*census_stage_*.ini"))
        if len(all_census) == 0:
            logging.error("No census stage found in the specified repository.")
            sys.exit()
        else:
            logging.info("%s census files have been detected.", str(len(all_census)))
        # Check that all census have the same stage
        all_census_stage = set(self.extract_census_stage(all_census))
        if len(all_census_stage) > 1:
            logging.error(
                "More than one type of census stage file have been detected: %s",
                ",".join(map(str, all_census_stage)),
            )
            logging.error("Potential mix between mapping output and profile output.")
            sys.exit()
        else:
            census_stage = list(all_census_stage)[0]
            assert census_stage in [1, 2]
            logging.info("Exactly one type of census found: %s", str(census_stage))

        # Create the dict: path -> ConfigParser
        all_census_dict = {
            my_census.parent: self.read_ini(my_census) for my_census in all_census
        }
        # Check there is exactly one census file per subdirectory
        try:
            assert len(all_census_dict) == len(all_census)
        except AssertionError:
            logging.error("There are several census files in the same directory. Exit.")
            sys.exit()
        # Create the dict: sample_name -> path
        all_samples_dict = {
            my_config["sample_info"]["sample_name"]: my_path
            for my_path, my_config in all_census_dict.items()
        }
        # Check there is exactly one sample per census file
        try:
            assert len(all_samples_dict) == len(all_census)
        except AssertionError:
            logging.error("Several census files refer to the same sample. Exit.")
            sys.exit()
        # Check that the ConfigParser match for mapping parameters
        logging.info("Checking that census parameters match...")
        if census_stage == 1:
            param_to_check = {"mapping": ["reference_name", "mapping_options"]}
        elif census_stage == 2:
            param_to_check = {
                "mapping": ["reference_name", "mapping_options", "database_type"],
                "profiling_parameters": [],
            }
        all_information = {
            my_path: self.extract_ini_info(my_config, param_to_check)
            for my_path, my_config in all_census_dict.items()
        }
        self.compare_section_info(all_information)
        # Save database_type for later use
        try:
            database_type_all = list(
                set(
                    [
                        my_info["database_type"]
                        for my_info in list(all_information.values())
                    ]
                )
            )
            assert len(database_type_all) == 1
            database_type = database_type_all[0]
        except:
            # Force to taxo if no information (e.g., census_stage_1) or no consensus
            database_type = "taxo"

        # Merge ini information
        logging.info("Merging ini information...")
        if census_stage == 1:
            param_to_save = {
                "mapping": [
                    "total_read_count",
                    "mapped_read_count",
                    "overall_alignment_rate",
                ]
            }
        elif census_stage == 2:
            param_to_save = {
                "mapping": [
                    "total_read_count",
                    "mapped_read_count",
                    "overall_alignment_rate",
                ],
                "profiling_stats": [],
            }
        all_information_to_save = {
            my_sample: self.extract_ini_info(all_census_dict[my_path], param_to_save)
            for my_sample, my_path in all_samples_dict.items()
        }
        all_information_to_save_df = (
            pd.DataFrame.from_dict(all_information_to_save)
            .transpose()
            .reset_index()
            .rename(columns={"index": "sample"})
        )
        output_name = self.output / f"{self.prefix}_report.tsv"
        all_information_to_save_df.to_csv(output_name, sep="\t", index=False)
        logging.info("Ini information saved as %s", output_name)

        # Merge the mapping output
        list_pattern_to_merge = {}
        if census_stage == 1:
            if not self.fast:
                list_pattern_to_merge = {".tsv": [self.meteor.gene_column]}
        elif census_stage == 2:
            list_pattern_to_merge["_msp.tsv"] = [self.meteor.msp_column]
            if not self.fast:
                list_pattern_to_merge["_genes.tsv"] = [self.meteor.gene_column]
            if database_type == "complete":
                list_pattern_to_merge["_modules_completeness.tsv"] = [
                    self.meteor.msp_column,
                    self.meteor.module_column,
                ]
                list_pattern_to_merge["_modules.tsv"] = [self.meteor.module_column]
                list_pattern_to_merge["_mustard.tsv"] = [self.meteor.ko_column]
        for my_pattern in list_pattern_to_merge:
            logging.info("Fetching output files with pattern %s", my_pattern)
            files_to_merge = self.find_files_to_merge(
                input_dir=all_samples_dict, pattern=my_pattern
            )
            logging.info(
                "There was %s files that correspond to the pattern.",
                len(files_to_merge),
            )
            logging.info("Merge the data.")
            merged_df = self.merge_df(files_to_merge, list_pattern_to_merge[my_pattern])
            logging.info("Save the data.")
            output_name = self.output / f"{self.prefix}{my_pattern}"
            merged_df.to_csv(output_name, sep="\t", index=False)
            logging.info("Data saved as %s", output_name)
