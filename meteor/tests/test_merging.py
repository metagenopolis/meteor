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

"""Test merging main object"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..merging import Merging
from pathlib import Path
import pytest
import pandas as pd
from configparser import ConfigParser


@pytest.fixture
def merging_profiles(datadir: Path, tmp_path: Path) -> Merging:
    meteor = Component
    meteor.profile_dir = datadir / "profiles"
    return Merging(meteor=meteor, output=tmp_path, prefix="my_test", fast=False)


@pytest.fixture
def merging_mapping(datadir: Path, tmp_path: Path) -> Merging:
    meteor = Component
    meteor.profile_dir = datadir / "mapping"
    return Merging(meteor=meteor, output=tmp_path, prefix="my_test", fast=False)


@pytest.fixture
def merging_fast(datadir: Path, tmp_path: Path) -> Merging:
    meteor = Component
    meteor.profile_dir = datadir / "profiles"
    return Merging(meteor=meteor, output=tmp_path, prefix="my_test", fast=True)


def test_extract_census_stage_1(merging_mapping: Merging) -> None:
    all_census = list(
        Path(merging_mapping.meteor.profile_dir).glob("**/*census_stage_*.ini")
    )
    all_census_stages = merging_mapping.extract_census_stage(all_census)
    assert all_census_stages == [1, 1, 1]


def test_extract_census_stage_2(merging_profiles: Merging) -> None:
    all_census = list(
        Path(merging_profiles.meteor.profile_dir).glob("**/*census_stage_*.ini")
    )
    all_census_stages = merging_profiles.extract_census_stage(all_census)
    assert all_census_stages == [2, 2, 2]


def test_find_files_to_merge(merging_profiles: Merging) -> None:
    path_dict = {
        "sample1": merging_profiles.meteor.profile_dir / "sample1",
        "sample2": merging_profiles.meteor.profile_dir / "sample2",
        "sample3": merging_profiles.meteor.profile_dir / "sample3",
    }
    list_files = merging_profiles.find_files_to_merge(path_dict, "_genes.tsv")
    assert list_files == {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_genes.tsv",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_genes.tsv",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_genes.tsv",
    }


def test_extract_ini_info(merging_profiles: Merging) -> None:
    config = ConfigParser()
    input_ini = (
        merging_profiles.meteor.profile_dir / "sample1" / "sample1_census_stage_2.ini"
    )
    with open(input_ini, "rt", encoding="UTF-8") as ini:
        config.read_file(ini)
    info = merging_profiles.extract_ini_info(
        config,
        param_dict={
            "profiling_parameters": ["msp_filter", "modules_def"],
            "mapping_file": [],
        },
    )
    assert info == {
        "msp_filter": "0.1",
        "modules_def": "GMM_definition.tsv",
        "bowtie_file_1": "sample1.sam",
        "mapping_file_format": "sam",
    }


def test_compare(merging_profiles: Merging) -> None:
    # Fetch all census ini files
    all_census = list(
        Path(merging_profiles.meteor.profile_dir).glob("**/*census_stage_*.ini")
    )
    # Create the dict: path -> ConfigParser
    all_census_dict = {
        my_census.parent: merging_profiles.read_ini(my_census)
        for my_census in all_census
    }
    # Define parameters that will be checked
    param_to_check = {
        "mapping": ["reference_name", "mapping_options", "database_type"],
        "profiling_parameters": [],
    }
    # Retrieve information about parameters
    all_information = {
        my_path: merging_profiles.extract_ini_info(my_config, param_to_check)
        for my_path, my_config in all_census_dict.items()
    }
    # Compare
    nb_inconsistencies = merging_profiles.compare_section_info(all_information)
    assert nb_inconsistencies == 1


def test_merge_df(merging_profiles: Merging, datadir: Path):
    # Test gene merging (1 key, same row numbers)
    files_to_merge = {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_genes.tsv",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_genes.tsv",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_genes.tsv",
    }
    merged_df = merging_profiles.merge_df(files_to_merge, key_merging=["gene_id"])
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_genes.tsv"
    )
    assert merged_df.equals(expected_output)
    # Test module merging (1 key, different row numbers)
    files_to_merge = {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_modules.tsv",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_modules.tsv",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_modules.tsv",
    }
    merged_df = merging_profiles.merge_df(files_to_merge, key_merging=["mod_id"])
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_modules.tsv"
    )
    assert merged_df.equals(expected_output)
    # Test module completeness merging (2 keys, different row numbers)
    files_to_merge = {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_modules_completeness.tsv",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_modules_completeness.tsv",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_modules_completeness.tsv",
    }
    merged_df = merging_profiles.merge_df(
        files_to_merge, key_merging=["msp_name", "mod_id"]
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_modules_completeness.tsv"
    )
    assert merged_df.equals(expected_output)


def test_execute1(merging_profiles: Merging, datadir: Path) -> None:
    merging_profiles.execute()

    # Check report
    real_output = merging_profiles.output / "my_test_report.tsv"
    assert real_output.exists()
    expected_output = (
        datadir / "expected_output" / "test_project_census_stage_2_report.tsv"
    )
    assert pd.read_table(real_output).equals(pd.read_table(expected_output))

    # Check existence and content of all files
    list_files = [
        "genes.tsv",
        "msp.tsv",
        "mustard.tsv",
        "modules.tsv",
        "modules_completeness.tsv",
    ]
    for my_file in list_files:
        real_output = merging_profiles.output / f"my_test_{my_file}"
        expected_output = datadir / "expected_output" / f"test_project_{my_file}"
        assert real_output.exists()
        real_output_df = pd.read_table(real_output)
        expected_output_df = pd.read_table(expected_output)
        real_output_df[["sample1", "sample2", "sample3"]] = real_output_df[
            ["sample1", "sample2", "sample3"]
        ].round(10)
        expected_output_df[["sample1", "sample2", "sample3"]] = expected_output_df[
            ["sample1", "sample2", "sample3"]
        ].round(10)
        assert real_output_df.equals(expected_output_df)


def test_execute2(merging_mapping: Merging, datadir: Path) -> None:
    merging_mapping.execute()

    # Check report
    real_output = merging_mapping.output / "my_test_report.tsv"
    assert real_output.exists()
    expected_output = (
        datadir / "expected_output" / "test_project_census_stage_1_report.tsv"
    )
    assert pd.read_table(real_output).equals(pd.read_table(expected_output))

    # Check existence and content of raw gene table
    real_output = merging_mapping.output / f"my_test.tsv"
    expected_output = datadir / "expected_output" / "test_project.tsv"
    assert real_output.exists()
    real_output_df = pd.read_table(real_output)
    expected_output_df = pd.read_table(expected_output)
    real_output_df[["sample1", "sample2", "sample3"]] = real_output_df[
        ["sample1", "sample2", "sample3"]
    ].round(10)
    expected_output_df[["sample1", "sample2", "sample3"]] = expected_output_df[
        ["sample1", "sample2", "sample3"]
    ].round(10)
    assert real_output_df.equals(expected_output_df)


def test_execute3(merging_fast: Merging, datadir: Path) -> None:
    merging_fast.execute()

    # Check report
    real_output = merging_fast.output / "my_test_report.tsv"
    assert real_output.exists()
    expected_output = (
        datadir / "expected_output" / "test_project_census_stage_2_report.tsv"
    )
    assert pd.read_table(real_output).equals(pd.read_table(expected_output))

    # Check existence and content of all files
    list_files = [
        "genes.tsv",
        "msp.tsv",
        "mustard.tsv",
        "modules.tsv",
        "modules_completeness.tsv",
    ]
    for my_file in list_files:
        real_output = merging_fast.output / f"my_test_{my_file}"
        expected_output = datadir / "expected_output" / f"test_project_{my_file}"
        if my_file == "genes.tsv":
            assert not real_output.exists()
        else:
            assert real_output.exists()
            real_output_df = pd.read_table(real_output)
            expected_output_df = pd.read_table(expected_output)
            real_output_df[["sample1", "sample2", "sample3"]] = real_output_df[
                ["sample1", "sample2", "sample3"]
            ].round(10)
            expected_output_df[["sample1", "sample2", "sample3"]] = expected_output_df[
                ["sample1", "sample2", "sample3"]
            ].round(10)
            assert real_output_df.equals(expected_output_df)
