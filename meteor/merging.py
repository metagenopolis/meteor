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

import numpy as np
import importlib.resources
import pandas as pd
import logging
import sys
import os
import fnmatch
from meteor.session import Session, Component
from meteor.qcreport import QCReport
from dataclasses import dataclass, field
from pathlib import Path
from biom.table import Table  # type: ignore
from typing import ClassVar
from functools import partial
from concurrent.futures import ThreadPoolExecutor


@dataclass
class Merging(Session):
    """
    Merge profiled outputs.
    """

    DEFAULT_MIN_MSP_ABUNDANCE: ClassVar[float] = 0.0
    DEFAULT_MIN_MSP_OCCURRENCE: ClassVar[int] = 1
    MPA_TAXONOMIC_LEVELS: ClassVar[list[str]] = [
        "a",
        "k",
        "p",
        "c",
        "o",
        "f",
        "g",
        "s",
        "t",
    ]
    DEFAULT_MPA_TAXONOMIC_LEVEL: ClassVar[str] = "a"
    DEFAULT_PREFIX: ClassVar[str] = "output"
    # Tables that can be merged: {table: {file pattern: key columns}}
    TABLES: ClassVar[dict[str, dict[str, list[str]]]] = {
        "msp": {"msp": ["msp_name"]},
        "genes": {"genes": ["gene_id"], "raw": ["gene_id"]},
        "modules_completeness": {"modules_completeness": ["msp_name", "mod_id"]},
        "modules": {"modules": ["mod_id"]},
        **{
            annotation: {
                f"{annotation}_as_genes_sum": ["annotation"],
                f"{annotation}_as_msp_sum": ["annotation"],
            }
            for annotation in ["mustard", "resfinder", "resfinderfg", "dbcan", "kegg"]
        },
    }
    # Tables of complete catalogues only
    FUNCTIONAL_TABLES: ClassVar[list[str]] = [
        table for table in TABLES if table not in ("msp", "genes")
    ]
    # Merged unless the tables are given (genes with -g)
    DEFAULT_TABLES: ClassVar[list[str]] = [table for table in TABLES if table != "genes"]

    meteor: type[Component]
    prefix: str
    min_msp_abundance: float
    min_msp_occurrence: int
    remove_sample_with_no_msp: bool
    output_mpa: bool
    mpa_taxonomic_level: str | None
    output_biom: bool
    output_gene_matrix: bool
    tables: list[str] | None = None
    sample_file: Path | None = None
    ranks: dict[str, str] = field(
        default_factory=lambda: {
            "superkingdom": "k__",
            "phylum": "p__",
            "class": "c__",
            "order": "o__",
            "family": "f__",
            "genus": "g__",
            "species": "s__",
            "msp_name": "t__",
        }
    )

    def __post_init__(self) -> None:
        self.meteor.merging_dir.mkdir(exist_ok=True, parents=True)

    # def extract_census_stage(self, census_list: List[Path]) -> List[int | None]:
    #     """For a list of census files formatted as
    #     sample1_census_stage_1.json or sample2_census_stage_2.json,
    #     extract the stage number (here [1, 2]).

    #     :param census_list: list of census_stage.json files as Path.
    #     """
    #     motif = re.compile(r"census_stage_(\d+)\.json")
    #     list_match = [motif.search(my_file.name) for my_file in census_list]
    #     return [int(my_match.group(1)) if my_match else None for my_match in list_match]

    def patterns_to_merge(self, database_type: str) -> dict[str, list[str]]:
        """Files to merge, from the tables asked (all but the genes by default,
        the genes with -g); functional tables only for complete catalogues.

        :param database_type: type of the catalogue (complete or taxo)
        :return: {file pattern: key columns}
        """
        selected = set(self.DEFAULT_TABLES if self.tables is None else self.tables)
        if self.output_gene_matrix:
            selected.add("genes")
        if database_type != "complete":
            functional = sorted(selected.intersection(self.FUNCTIONAL_TABLES))
            if functional and self.tables is not None:
                logging.error(
                    "The %s tables can only be merged with a complete catalogue.",
                    ", ".join(functional),
                )
                sys.exit(1)
            selected.difference_update(self.FUNCTIONAL_TABLES)
        if (self.output_biom or self.output_mpa) and "msp" not in selected:
            logging.warning("The species table is not merged: no biom or mpa file.")
        return {
            pattern: keys
            for table, patterns in self.TABLES.items()
            if table in selected
            for pattern, keys in patterns.items()
        }

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
        len_list = list({len(value) for value in list(dict_to_merge.values())})
        assert len(len_list) == 1
        assert len_list[0] == 1
        files_to_merge = {
            my_sample: my_list[0] for my_sample, my_list in dict_to_merge.items()
        }
        return files_to_merge

    def extract_json_info(
        self, config: dict, param_dict: dict[str, list[str]]
    ) -> dict[str, str]:
        """Get information from the required section and fields

        :param config: A dict from provided json file.
        :param param_dict: A dictionnary in the form of {section1: [fields1, fields2], section2: []}
        giving the information that should be extracted from the config.
        """
        # Check that sections are present
        try:
            assert all(my_section in config for my_section in param_dict.keys())
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
                my_field in config[my_section]
                for my_section in param_dict
                for my_field in param_dict[my_section]
            )
        except AssertionError:
            logging.error("Missing required fields in census json file.")
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
        logging.info("Use the json file associated to %s as a reference.", ref_config)
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
                "There was %d files that do not match the reference json file",
                len(problematic_config),
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

        :return: merged data frame
        """
        # Load the data frames
        list_df = [
            pd.read_table(my_path, compression="xz")
            .rename(columns={"value": my_sample})[key_merging + [my_sample]]
            .assign(
                **{
                    my_sample: lambda df: df[my_sample].astype(
                        pd.SparseDtype("float", 0.0)
                    )
                }
            )
            for my_sample, my_path in dict_path.items()
        ]
        for df in list_df:
            df.set_index(key_merging, inplace=True)

        # Now you can concatenate using the aligned index
        merged_df = pd.concat(list_df, join="outer", axis=1).reset_index()
        return merged_df

    def _write_annotations(
        self, my_pattern: str, filtered_df: pd.DataFrame, output_name: Path
    ) -> None:
        """Annotation tables of the kept modules / mustard / kegg rows

        :param my_pattern: merged profile pattern
        :param filtered_df: kept rows (key columns are enough)
        :param output_name: output prefix
        """
        if my_pattern == "modules":
            module_path = Path(
                str(
                    importlib.resources.files("meteor")
                    / "data/modules_definition.feather"
                )
            )
            annotation = self.load_data(module_path)
            annotation = annotation[annotation["id"].isin(filtered_df["mod_id"])]
            annotation[
                [
                    "id",
                    "type",
                    "name",
                ]
            ].to_csv(
                Path(f"{output_name}_definition.tsv"),
                sep="\t",
                index=False,
            )
        elif my_pattern == "mustard_as_genes_sum":
            mustard_annotation_path = Path(
                str(
                    importlib.resources.files("meteor")
                    / "data/category_pcm.feather"
                )
            )
            mustard_annotation = self.load_data(mustard_annotation_path)
            mustard_annotation = mustard_annotation[
                mustard_annotation["ARD"].isin(filtered_df["annotation"])
            ]
            mustard_annotation[
                [
                    "ARD",
                    "Antimicrobial"]
            ].to_csv(
                Path(f"{output_name}_antimicrobial.tsv"),
                sep="\t",
                index=False,
            )
        elif my_pattern == "kegg_as_genes_sum":
            kegg_annotation_path = Path(
                str(
                    importlib.resources.files("meteor")
                    / "data/ko_list_v116.feather"
                )
            )
            kegg_annotation = self.load_data(kegg_annotation_path)
            kegg_annotation = kegg_annotation[
                kegg_annotation["knum"].isin(filtered_df["annotation"])
            ]
            kegg_annotation[
                [
                    "knum",
                    "definition"]
            ].to_csv(
                Path(f"{output_name}_description.tsv"),
                sep="\t",
                index=False,
            )

    # rows formatted and written per block by _merge_filter_write
    WRITE_BLOCK_ROWS: ClassVar[int] = 4096

    def mpa_table(
        self, msp_df: pd.DataFrame, taxonomy: pd.DataFrame
    ) -> pd.DataFrame:
        """Species abundance table in the style of MetaPhlAn (merged profiles).

        One row per clade, from k__ (superkingdom) to t__ (MSP), named by its
        lineage (k__Bacteria|p__Bacillota|...|s__Escherichia_coli|t__msp_0005);
        the abundance of a clade is the sum of its MSPs, in percent of the
        sample (each rank sums to 100). Clades absent from every sample are not
        written. Ranks missing from the taxonomy are named "unclassified".

        :param msp_df: merged MSP table (msp_name, then one column per sample)
        :param taxonomy: catalogue taxonomy (msp_name and the ranks of self.ranks)
        :return: table indexed by clade_name, one column per sample
        """
        samples = msp_df.columns.drop("msp_name")
        abundance = np.nan_to_num(msp_df[samples].to_numpy(dtype=float))
        totals = abundance.sum(axis=0)
        percent = pd.DataFrame(
            100 * abundance / np.where(totals > 0, totals, 1.0),
            index=msp_df["msp_name"].to_numpy(),
            columns=samples,
        )
        taxonomic_ranks = [rank for rank in self.ranks if rank != "msp_name"]
        lineage = (
            taxonomy.drop_duplicates("msp_name")
            .set_index("msp_name")[taxonomic_ranks]
            .reindex(percent.index)
        )
        not_annotated = lineage.isna().all(axis=1)
        if not_annotated.any():
            logging.warning(
                "%d MSPs are not in the taxonomy of the catalogue, they are "
                "reported as unclassified in the mpa table.",
                not_annotated.sum(),
            )
        lineage = lineage.fillna("").astype(str)
        # MetaPhlAn clade names: no space, no "|"
        lineage = lineage.apply(
            lambda names: names.str.strip()
            .str.replace(r"[\s|]+", "_", regex=True)
            .replace("", "unclassified")
        )
        lineage["msp_name"] = lineage.index
        levels = (
            self.MPA_TAXONOMIC_LEVELS[1:]
            if self.mpa_taxonomic_level in (None, "a")
            else [self.mpa_taxonomic_level]
        )
        clades = []
        clade_name = pd.Series("", index=lineage.index)
        for rank, prefix in self.ranks.items():
            clade_name = clade_name + ("|" if rank != "superkingdom" else "") + prefix
            clade_name = clade_name + lineage[rank]
            if prefix[0] in levels:
                clades.append(percent.groupby(clade_name.to_numpy(), sort=False).sum())
        table = pd.concat(clades)
        table = table.loc[table.sum(axis=1) > 0]
        # MetaPhlAn order: by rank, then by name
        table = table.loc[sorted(table.index, key=lambda name: (name.count("|"), name))]
        table.index.name = "clade_name"
        return table

    def write_mpa(
        self, table: pd.DataFrame, ref_json: dict, output_file: Path
    ) -> None:
        """Write the mpa table with a first line naming meteor, catalogue and
        taxonomy versions, as merge_metaphlan_tables.py writes the MetaPhlAn
        database version.

        :param table: mpa_table output
        :param ref_json: catalogue reference json
        :param output_file: output tsv file
        """
        header = f"#meteor_v{self.meteor.version}"
        reference_name = ref_json.get("reference_info", {}).get("reference_name")
        if reference_name:
            header += f"_{reference_name}"
        taxonomy_version = ref_json["annotation"]["taxonomy"].get("version")
        if taxonomy_version:
            header += f"_gtdb_r{taxonomy_version}"
        with open(output_file, "w", encoding="utf-8") as out:
            out.write(header + "\n")
            table.to_csv(out, sep="\t")

    def _merge_filter_write(
        self, dict_path: dict[str, Path], key_merging: list[str], output_file: Path
    ) -> pd.DataFrame | None:
        """merge_df + abundance/occurrence filters + sample filter + to_csv without
        a (sparse) DataFrame holding every sample.

        Writes the same table as the DataFrame path of execute and returns the key
        columns of the kept rows, or None when the table must be built with the
        DataFrame path (duplicated keys, no row kept, field needing quotes).

        :param dict_path: {sample name: profile file}, in output column order
        :param key_merging: key columns
        :param output_file: output tsv file
        """
        def read(my_path: Path) -> tuple[pd.Index, np.ndarray]:
            df = pd.read_table(my_path, compression="xz")
            if len(key_merging) > 1:
                index = pd.MultiIndex.from_frame(df[key_merging])
            else:
                index = pd.Index(df[key_merging[0]], name=key_merging[0])
            return index, df["value"].to_numpy(dtype=float)

        # Files are independent: read them in parallel (decompression and
        # parsing release the GIL); results are kept in sample order
        try:
            n_cpus = len(os.sched_getaffinity(0))
        except AttributeError:
            n_cpus = os.cpu_count() or 1
        with ThreadPoolExecutor(max_workers=max(1, min(8, n_cpus))) as executor:
            loaded = list(executor.map(read, dict_path.values()))
        indexes = [index for index, _ in loaded]
        values = [vals for _, vals in loaded]
        del loaded
        # Same row order as pd.concat(..., join="outer", axis=1)
        try:
            # what pd.concat uses to combine the row indexes
            from pandas.core.indexes.api import _get_combined_index

            union = _get_combined_index(indexes, intersect=False, sort=False)
            # pandas >= 3.1 returns (index, whether all indexes are equal)
            if isinstance(union, tuple):
                union = union[0]
            if not isinstance(union, pd.Index):
                raise TypeError(f"unexpected combined index: {type(union)}")
        except Exception:  # pylint: disable=broad-except
            try:
                union = pd.concat(
                    [pd.DataFrame(index=index) for index in indexes], axis=1, join="outer"
                ).index
            except Exception:  # pylint: disable=broad-except
                return None
        if not union.is_unique:
            return None
        # Column-major matrix, as DataFrame.to_numpy() of the merged sparse columns
        matrix = np.full((len(union), len(indexes)), np.nan, order="F")
        for j, (index, vals) in enumerate(zip(indexes, values)):
            positions = union.get_indexer(index)
            if np.bincount(positions, minlength=len(union)).max(initial=0) > 1:
                # duplicated keys in a sample: the DataFrame path raises as pandas does
                return None
            matrix[positions, j] = vals
        del indexes, values
        # Row filters, computed on blocks of rows to avoid full-size temporaries
        # (each row is reduced the same way as on the whole matrix)
        row_sums = np.empty(len(union))
        occurrence = np.empty(len(union), dtype=np.int64)
        for start in range(0, len(union), 65536):
            block = matrix[start : start + 65536]
            row_sums[start : start + 65536] = np.nansum(block, axis=1)
            occurrence[start : start + 65536] = np.nansum(
                (block != 0) & ~np.isnan(block), axis=1
            )
        kept_rows = np.flatnonzero(
            (row_sums >= self.min_msp_abundance) & (occurrence >= self.min_msp_occurrence)
        )
        if len(kept_rows) == 0:
            return None
        samples = list(dict_path)
        kept_cols = np.arange(len(samples))
        if self.remove_sample_with_no_msp:
            # key columns are kept (non-empty), samples with a non-zero sum are kept
            col_sums = np.array(
                [np.nansum(matrix[kept_rows, j]) for j in range(len(samples))]
            )
            kept_cols = np.flatnonzero(col_sums != 0)
        keys = union.to_frame(index=False) if isinstance(union, pd.MultiIndex) \
            else pd.DataFrame({key_merging[0]: union})
        keys = keys.iloc[kept_rows].reset_index(drop=True)
        if keys.isna().to_numpy().any():
            return None
        header = key_merging + [samples[j] for j in kept_cols]
        key_str = keys.astype(str).to_numpy()
        special = ("\t", '"', "\n", "\r")
        if any(c in str(field) for field in header for c in special) or any(
            np.char.find(key_str.astype(str), c).max() >= 0 for c in special
        ):
            return None
        with open(output_file, "w", encoding="utf-8") as out:
            out.write("\t".join(map(str, header)) + "\n")
            for start in range(0, len(kept_rows), self.WRITE_BLOCK_ROWS):
                rows = kept_rows[start : start + self.WRITE_BLOCK_ROWS]
                block = matrix[np.ix_(rows, kept_cols)]
                # values written as DataFrame.to_csv does: str() of the float,
                # missing values as empty fields
                uniques, inverse = np.unique(block, return_inverse=True)
                text = uniques.astype(str).astype(object)
                text[np.isnan(uniques)] = ""
                cells = text[inverse.reshape(block.shape)]
                lines = [
                    "\t".join(k) + "\t" + "\t".join(c) if len(c) else "\t".join(k)
                    for k, c in zip(key_str[start : start + len(rows)].tolist(), cells.tolist())
                ]
                out.write("\n".join(lines) + "\n")
        return keys

    def execute(self) -> None:
        "Merge all files generated by either profiler."
        # Load reference data
        ref_json = self.get_reference_info(self.meteor.ref_dir)
        database_type = ref_json["reference_info"]["database_type"]
        #Fetch all census ini files
        all_census = list(Path(self.meteor.profile_dir).glob("**/*census_stage_2.json"))
        all_census = [
            Path(root, filename)
            for root, _, files in os.walk(self.meteor.profile_dir, followlinks=True)
            for filename in fnmatch.filter(files, "*census_stage_2.json")
        ]
        if len(all_census) == 0:
            logging.error(
                "No *_census_stage_2.json files found in %s",
                self.meteor.profile_dir,
            )
            sys.exit(1)
        else:
            logging.info("%d census files have been detected.", len(all_census))
        # Create the dict: path -> Dict
        all_census_dict = {
            my_census.parent: self.read_json(my_census) for my_census in all_census
        }
        # # The meteor version is not a constraint to matching
        # print(all_census_dict)
        # for my_census in all_census_dict:
        #     del all_census_dict[my_census]["meteor_version"]
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
        # Samples of the sample list (-l), in its order
        all_samples_dict = self.select_samples(all_samples_dict, self.sample_file)
        all_census_dict = {
            my_path: all_census_dict[my_path] for my_path in all_samples_dict.values()
        }
        # Check that the json match for mapping parameters
        logging.info("Checking that census parameters match...")
        param_to_check = {
            "mapping": [
                "reference_name",
                "trim",
                "alignment_number",
                "mapping_type",
                "database_type",
            ],
            "counting": [
                "identity_threshold",
            ],
            "profiling_parameters": [""],
        }
        all_information = {
            my_path: self.extract_json_info(my_config, param_to_check)
            for my_path, my_config in all_census_dict.items()
        }
        self.compare_section_info(all_information)

        # Merge json information
        logging.info("Merging json information...")
        # Get all values from all fields from all sections from all json files
        all_information_to_save = {
            key.name: {
                option: value
                for section in config_parser
                if section != "meteor_version"
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
        output_name = self.meteor.merging_dir / f"{self.prefix}_report.tsv"
        all_information_to_save_df.to_csv(output_name, sep="\t", index=False)
        logging.info("Json information saved as %s", output_name)

        # Merge the mapping output
        list_pattern_to_merge = self.patterns_to_merge(database_type)
        msp_merged, taxonomy = None, None
        for my_pattern, value in list_pattern_to_merge.items():
            logging.info(
                "Fetching output files with pattern %s", my_pattern + ".tsv.xz"
            )
            files_to_merge = self.find_files_to_merge(
                input_dir=all_samples_dict, pattern=f"{my_pattern}.tsv.xz"
            )
            logging.info(
                "There was %s files that correspond to the pattern.",
                len(files_to_merge),
            )
            output_name = self.meteor.merging_dir / f"{self.prefix}_{my_pattern}"
            if my_pattern != "msp":
                # Large tables: merged, filtered and written without the sparse
                # DataFrame (same output file)
                logging.info("Merge and save the data.")
                kept_keys = self._merge_filter_write(
                    files_to_merge, value, output_name.with_suffix(".tsv")
                )
                if kept_keys is not None:
                    self._write_annotations(my_pattern, kept_keys, output_name)
                    logging.info("Data saved as %s", output_name)
                    continue
            logging.info("Merge the data.")
            merged_df = self.merge_df(files_to_merge, value)
            logging.info("Save the data.")
            # Calculate the occurrence of non-zero values across the specified columns only
            # occurrence = (merged_df[list_pattern_to_merge[my_pattern]] != 0).sum(axis=1)
            # Calculate the sum of each row in the DataFrame
            numeric_df = merged_df.drop(columns=value).to_numpy()
            # row_sums = numeric_df.sum(axis=1)
            row_sums = np.nansum(numeric_df, axis=1)

            # Calculate the occurrence (number of non-zero values) in each row
            # occurrence = (numeric_df != 0).sum(axis=1)
            occurrence = np.nansum((numeric_df != 0) & ~np.isnan(numeric_df), axis=1)
            # Apply both abundance and occurrence filters directly
            # print(merged_df[list_pattern_to_merge[my_pattern]])
            # Filter the DataFrame based on the calculated row sums and occurrence
            filtered_df = merged_df.loc[
                (row_sums >= self.min_msp_abundance)
                & (occurrence >= self.min_msp_occurrence),
                :,
            ]
            if self.remove_sample_with_no_msp:
                # Filter the DataFrame to keep only columns with a non-zero sum
                filtered_df = filtered_df.loc[:, (filtered_df.sum(axis=0) != 0)]
            if my_pattern == "msp":
                annotation = self.load_data(
                    self.meteor.ref_dir
                    / ref_json["reference_file"]["database_dir"]
                    / ref_json["annotation"]["taxonomy"]["filename"]
                )
                # species table before filters, for the report
                msp_merged, taxonomy = merged_df, annotation
                annotation = annotation[
                    annotation["msp_name"].isin(filtered_df["msp_name"])
                ]
                annotation[
                    [
                        "msp_name",
                        "superkingdom",
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ]
                ].to_csv(
                    Path(f"{output_name}_taxonomy.tsv"),
                    sep="\t",
                    index=False,
                )
                if self.output_mpa:
                    mpa_output = Path(f"{output_name}_mpa.tsv")
                    self.write_mpa(
                        self.mpa_table(filtered_df, annotation), ref_json, mpa_output
                    )
                    logging.info("MetaPhlAn-style table saved as %s", mpa_output)
                if self.output_biom:
                    sample_ids = filtered_df.columns.drop("msp_name").tolist()
                    observ_ids = filtered_df["msp_name"].tolist()
                    # Apply the prefixes to each taxonomic rank
                    for rank, prefix in self.ranks.items():
                        annotation[rank] = annotation[rank].apply(
                            partial(lambda prefix, x: f"{prefix}{x}", prefix)
                        )
                    observ_metadata = [
                        {"taxonomy": row.iloc[1:].tolist()}
                        for _, row in annotation.iterrows()
                    ]
                    biom_table = Table(
                        data=filtered_df.drop(columns=["msp_name"]).to_numpy(),
                        observation_ids=observ_ids,
                        sample_ids=sample_ids,
                        observation_metadata=observ_metadata,
                        # table_id="Meteor Table",
                        type="Taxon table",
                    )
                    # Generate JSON representation of the BIOM table
                    biom_json = biom_table.to_json(generated_by="Meteor")
                    # Write the JSON to a file
                    with open(
                        output_name.with_suffix(".biom"), "wt", encoding="UTF-8"
                    ) as f:
                        f.write(biom_json)
                    # with h5py.File(output_name.with_suffix(".biom"), "w") as f:
                    #     table.to_hdf5(f, generated_by="Meteor", compress=True)
            else:
                self._write_annotations(my_pattern, filtered_df, output_name)
            filtered_df.to_csv(output_name.with_suffix(".tsv"), sep="\t", index=False)
            logging.info("Data saved as %s", output_name)
        # Mapped reads, MSP richness and composition per sample, with MultiQC files
        QCReport(self.meteor.merging_dir, self.prefix).execute(
            {name: all_census_dict[path] for name, path in all_samples_dict.items()},
            msp_merged,
            taxonomy,
        )
