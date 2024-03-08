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
from dataclasses import dataclass
from typing import Type, Dict, List
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
    output: Path
    prefix: str
    fast: bool = False

    def __post_init__(self) -> None:
        self.output.mkdir(exist_ok=True, parents=True)

    # def extract_census_stage(self, census_list: List[Path]) -> List[int | None]:
    #     """For a list of census files formatted as
    #     sample1_census_stage_1.json or sample2_census_stage_2.json,
    #     extract the stage number (here [1, 2]).

    #     :param census_list: list of census_stage.json files as Path.
    #     """
    #     motif = re.compile(r"census_stage_(\d+)\.json")
    #     list_match = [motif.search(my_file.name) for my_file in census_list]
    #     return [int(my_match.group(1)) if my_match else None for my_match in list_match]

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

    def extract_json_info(
        self, config: Dict, param_dict: Dict[str, List[str]]
    ) -> dict[str, str]:
        """Get information from the required section and fields

        :param config: A dict from provided json file.
        :param param_dict: A dictionnary in the form of {section1: [fields1, fields2], section2: []}
        giving the information that should be extracted from the config.
        """
        # Check that sections are present
        try:
            assert all([my_section in config for my_section in list(param_dict.keys())])
        except AssertionError:
            logging.error("Missing required section in census json file.")
            sys.exit(1)
        # If fields is empty, consider all fields of the section
        param_dict = {
            key: list(config[key]) if value == [""] else value
            for key, value in param_dict.items()
        }
        # Check that required fields are present
        try:
            assert all(
                [
                    my_field in config[my_section]
                    for my_section in param_dict
                    for my_field in param_dict[my_section]
                ]
            )
        except AssertionError:
            logging.error("Missing required fields in census ini file.")
            sys.exit(1)
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
        ref_config = sorted(list(all_sections_info.keys()))[0]
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
                "There was %s files that do not match the reference json file",
                str(len(problematic_config)),
            )
            for i in problematic_config:
                logging.warning(
                    "The file %s does not match the reference json file.", i
                )
        else:
            logging.info("All files have the same parameters.")
        return len(problematic_config)

    def merge_df(
        self, dict_path: dict[str, Path], key_merging: list[str]
    ) -> pd.DataFrame:
        """Merge data frames contained in a list according to a pivot key.

        :param dict_path: dictionnary of path of files to merge, where key is sample name
        and will be used as colnames.
        :param key_merging: list of keys that should be used as pivot for merging.
        """
        # Load the data frames
        list_df = [
            pd.read_table(my_path, compression="xz").rename(
                columns={"value": my_sample}
            )[key_merging + [my_sample]]
            for (my_sample, my_path) in dict_path.items()
        ]
        merged_df = reduce(
            lambda left, right: pd.merge(left, right, how="outer"),
            list_df,
        )
        return merged_df

    def execute(self) -> None:
        "Merge all files generated by either profiler."
        # Fetch all census ini files
        all_census = list(Path(self.meteor.profile_dir).glob("**/*census_stage_2.json"))
        if len(all_census) == 0:
            logging.error("No census stage 2 found in the specified repository.")
            sys.exit(1)
        else:
            logging.info("%s census files have been detected.", str(len(all_census)))
        # Create the dict: path -> Dict
        all_census_dict = {
            my_census.parent: self.read_json(my_census) for my_census in all_census
        }
        # Check there is exactly one census file per subdirectory
        try:
            assert len(all_census_dict) == len(all_census)
        except AssertionError:
            logging.error("There are several census files in the same directory.")
            sys.exit(1)
        # Create the dict: sample_name -> path
        all_samples_dict = {
            my_config["sample_info"]["sample_name"]: my_path
            for my_path, my_config in all_census_dict.items()
        }
        # Check there is exactly one sample per census file
        try:
            assert len(all_samples_dict) == len(all_census)
        except AssertionError:
            logging.error("Several census files refer to the same sample.")
            sys.exit(1)
        # Check that the json match for mapping parameters
        logging.info("Checking that census parameters match...")
        param_to_check = {
            "mapping": [
                "reference_name",
                "trim",
                "alignment_number",
                "mapping_type",
                "identity_threshold",
                "database_type",
            ],
            "profiling_parameters": [""],
        }
        all_information = {
            my_path: self.extract_json_info(my_config, param_to_check)
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
        except AssertionError:
            # Force to taxo in no consensus
            database_type = "taxo"

        # Merge ini information
        logging.info("Merging json information...")
        # Get all values from all fields from all sections from all json files
        all_information_to_save = {
            key.name: {
                option: value
                for section in config_parser
                for option, value in config_parser[section].items()
            }
            for key, config_parser in all_census_dict.items()
        }

        all_information_to_save_df = (
            pd.DataFrame.from_dict(all_information_to_save)
            .transpose()
            .reset_index()
            .rename(columns={"index": "sample"})
        )
        output_name = self.output / f"{self.prefix}_report.tsv"
        all_information_to_save_df.to_csv(output_name, sep="\t", index=False)
        logging.info("Json information saved as %s", output_name)

        # Merge the mapping output
        list_pattern_to_merge = {}
        list_pattern_to_merge["_msp.tsv"] = ["msp_name"]
        if not self.fast:
            list_pattern_to_merge["_genes.tsv"] = ["gene_id"]
            list_pattern_to_merge["_raw.tsv"] = ["gene_id"]
        if database_type == "complete":
            list_pattern_to_merge["_modules_completeness.tsv"] = [
                "msp_name",
                "mod_id",
            ]
            list_pattern_to_merge["_modules.tsv"] = ["mod_id"]
            list_pattern_to_merge["_mustard_as_genes_sum.tsv"] = ["annotation"]
            list_pattern_to_merge["_dbcan_as_msp_sum.tsv"] = ["annotation"]
        for my_pattern in list_pattern_to_merge:
            logging.info("Fetching output files with pattern %s", my_pattern + ".xz")
            files_to_merge = self.find_files_to_merge(
                input_dir=all_samples_dict, pattern=my_pattern + ".xz"
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
