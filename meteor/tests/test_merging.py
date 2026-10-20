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
import json


@pytest.fixture
def merging_profiles(datadir: Path, tmp_path: Path) -> Merging:
    meteor = Component
    meteor.profile_dir = datadir / "profiles"
    meteor.merging_dir = tmp_path
    meteor.ref_dir = datadir / "ref_dir"
    return Merging(
        meteor=meteor,
        prefix="my_test",
        min_msp_abundance=0.0,
        min_msp_occurrence=0,
        remove_sample_with_no_msp=False,
        output_mpa=False,
        mpa_taxonomic_level=None,
        output_biom=False,
        output_gene_matrix=True,
    )


# @pytest.fixture
# def merging_mapping(datadir: Path, tmp_path: Path) -> Merging:
#     meteor = Component
#     meteor.profile_dir = datadir / "mapping"
#     return Merging(meteor=meteor, output=tmp_path, prefix="my_test", fast=False)


@pytest.fixture
def merging_fast(datadir: Path, tmp_path: Path) -> Merging:
    meteor = Component
    meteor.profile_dir = datadir / "profiles"
    meteor.merging_dir = tmp_path
    meteor.ref_dir = datadir / "ref_dir"
    return Merging(
        meteor=meteor,
        prefix="my_test",
        min_msp_abundance=0.0,
        min_msp_occurrence=0,
        remove_sample_with_no_msp=False,
        output_mpa=False,
        mpa_taxonomic_level=None,
        output_biom=False,
        output_gene_matrix=False,
    )


# def test_extract_census_stage_1(merging_mapping: Merging) -> None:
#     all_census = list(
#         Path(merging_mapping.meteor.profile_dir).glob("**/*census_stage_*.json")
#     )
#     all_census_stages = merging_mapping.extract_census_stage(all_census)
#     assert all_census_stages == [1, 1, 1]


# def test_extract_census_stage_2(merging_profiles: Merging) -> None:
#     all_census = list(
#         Path(merging_profiles.meteor.profile_dir).glob("**/*census_stage_*.json")
#     )
#     all_census_stages = merging_profiles.extract_census_stage(all_census)
#     assert all_census_stages == [2, 2, 2]


def test_find_files_to_merge(merging_profiles: Merging) -> None:
    path_dict = {
        "sample1": merging_profiles.meteor.profile_dir / "sample1",
        "sample2": merging_profiles.meteor.profile_dir / "sample2",
        "sample3": merging_profiles.meteor.profile_dir / "sample3",
    }
    list_files = merging_profiles.find_files_to_merge(path_dict, "_genes.tsv.xz")
    assert list_files == {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_genes.tsv.xz",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_genes.tsv.xz",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_genes.tsv.xz",
    }


def test_extract_json_info(merging_profiles: Merging) -> None:
    config = {}
    input_json = (
        merging_profiles.meteor.profile_dir / "sample1" / "sample1_census_stage_2.json"
    )
    with open(input_json, "rt", encoding="UTF-8") as json_data:
        config = json.load(json_data)
    info = merging_profiles.extract_json_info(
        config,
        param_dict={
            "profiling_parameters": ["msp_filter", "modules_def"],
            "mapping": ["mapping_file"],
        },
    )
    assert info == {
        "msp_filter": 0.1,
        "modules_def": "modules_definition.tsv",
        "mapping_file": "sample1.sam",
    }


def test_compare(merging_profiles: Merging) -> None:
    # Fetch all census ini files
    all_census = list(
        Path(merging_profiles.meteor.profile_dir).glob("**/*census_stage_*.json")
    )
    # Create the dict: path -> Dict
    all_census_dict = {
        my_census.parent: merging_profiles.read_json(my_census)
        for my_census in all_census
    }
    # Define parameters that will be checked
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
    # Retrieve information about parameters
    all_information = {
        my_path: merging_profiles.extract_json_info(my_config, param_to_check)
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
        / "sample1_genes.tsv.xz",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_genes.tsv.xz",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_genes.tsv.xz",
    }
    merged_df = merging_profiles.merge_df(files_to_merge, key_merging=["gene_id"])
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_genes.tsv"
    )
    for col in expected_output.select_dtypes(include=["float64"]).columns:
        expected_output[col] = expected_output[col].astype(pd.SparseDtype(float, 0.0))
    assert merged_df.equals(expected_output)
    # Test module merging (1 key, different row numbers)
    files_to_merge = {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_modules.tsv.xz",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_modules.tsv.xz",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_modules.tsv.xz",
    }
    merged_df = merging_profiles.merge_df(files_to_merge, key_merging=["mod_id"])
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_modules.tsv"
    )
    for col in expected_output.select_dtypes(include=["float64"]).columns:
        expected_output[col] = expected_output[col].astype(pd.SparseDtype(float, 0.0))
    assert merged_df.equals(expected_output)
    # Test module completeness merging (2 keys, different row numbers)
    files_to_merge = {
        "sample1": merging_profiles.meteor.profile_dir
        / "sample1"
        / "sample1_modules_completeness.tsv.xz",
        "sample2": merging_profiles.meteor.profile_dir
        / "sample2"
        / "sample2_modules_completeness.tsv.xz",
        "sample3": merging_profiles.meteor.profile_dir
        / "sample3"
        / "sample3_modules_completeness.tsv.xz",
    }
    merged_df = merging_profiles.merge_df(
        files_to_merge, key_merging=["msp_name", "mod_id"]
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "test_project_modules_completeness.tsv"
    )
    for col in expected_output.select_dtypes(include=["float64"]).columns:
        expected_output[col] = expected_output[col].astype(pd.SparseDtype(float, 0.0))
    merged_df = merged_df.sort_values(by=["msp_name", "mod_id"]).reset_index(drop=True)
    assert merged_df.equals(expected_output)


def test_execute1(merging_profiles: Merging, datadir: Path) -> None:
    merging_profiles.execute()

    # What is in the directory
    for path in Path(merging_profiles.meteor.merging_dir).iterdir():
        print(path)

    # Check report
    real_output = merging_profiles.meteor.merging_dir / "my_test_report.tsv"
    assert real_output.exists()
    expected_output = (
        datadir / "expected_output" / "test_project_census_stage_2_report.tsv"
    )
    real_output_df = pd.read_table(real_output)

    expected_output_df = pd.read_table(expected_output)
    real_output_df = (
        real_output_df.sort_values(by=["sample"])
        .reset_index(drop=True)
        .reindex(sorted(real_output_df.columns), axis=1)
    )
    expected_output_df = (
        expected_output_df.sort_values(by=["sample"])
        .reset_index(drop=True)
        .reindex(sorted(expected_output_df.columns), axis=1)
    )
    assert real_output_df.round(2).equals(expected_output_df.round(2))

    # Check existence and content of all files
    list_files = [
        "raw.tsv",
        "genes.tsv",
        "msp.tsv",
        "mustard_as_genes_sum.tsv",
        "dbcan_as_msp_sum.tsv",
        "modules.tsv",
        "modules_completeness.tsv",
    ]
    for my_file in list_files:
        real_output = merging_profiles.meteor.merging_dir / f"my_test_{my_file}"
        expected_output = datadir / "expected_output" / f"test_project_{my_file}"
        assert real_output.exists()
        real_output_df = pd.read_table(real_output).reindex(
            sorted(real_output_df.columns), axis=1
        )
        expected_output_df = pd.read_table(expected_output).reindex(
            sorted(expected_output_df.columns), axis=1
        )
        assert real_output_df.round(10).equals(expected_output_df.round(10))


def test_execute2(merging_fast: Merging, datadir: Path) -> None:
    merging_fast.execute()

    # Check report
    real_output = merging_fast.meteor.merging_dir / "my_test_report.tsv"
    assert real_output.exists()
    expected_output = (
        datadir / "expected_output" / "test_project_census_stage_2_report.tsv"
    )
    real_output_df = pd.read_table(real_output)
    expected_output_df = pd.read_table(expected_output)
    real_output_df = (
        real_output_df.sort_values(by=["sample"])
        .reset_index(drop=True)
        .reindex(sorted(real_output_df.columns), axis=1)
    )
    expected_output_df = (
        expected_output_df.sort_values(by=["sample"])
        .reset_index(drop=True)
        .reindex(sorted(expected_output_df.columns), axis=1)
    )
    assert real_output_df.round(2).equals(expected_output_df.round(2))

    # Check existence and content of all files
    list_files = [
        "raw.tsv",
        "genes.tsv",
        "msp.tsv",
        "mustard_as_genes_sum.tsv",
        "dbcan_as_msp_sum.tsv",
        "modules.tsv",
        "modules_completeness.tsv",
    ]
    for my_file in list_files:
        real_output = merging_fast.meteor.merging_dir / f"my_test_{my_file}"
        expected_output = datadir / "expected_output" / f"test_project_{my_file}"
        if my_file in ["genes.tsv", "raw.tsv"]:
            assert not real_output.exists()
        else:
            assert real_output.exists()
            real_output_df = pd.read_table(real_output).reindex(
                sorted(real_output_df.columns), axis=1
            )
            expected_output_df = pd.read_table(expected_output).reindex(
                sorted(expected_output_df.columns), axis=1
            )
            assert real_output_df.round(10).equals(expected_output_df.round(10))


@pytest.mark.parametrize("remove_samples", [False, True])
def test_merge_filter_write_matches_dataframe_path(
    merging_profiles: Merging, tmp_path: Path, remove_samples: bool
) -> None:
    """_merge_filter_write writes the same table as merge_df + filters + to_csv"""
    import numpy as np
    merging_profiles.remove_sample_with_no_msp = remove_samples
    merging_profiles.min_msp_occurrence = 1
    samples = {
        name: merging_profiles.meteor.profile_dir / name for name in ("sample1", "sample2", "sample3")
    }
    for pattern, keys in [("genes", ["gene_id"]), ("modules", ["mod_id"]),
                          ("modules_completeness", ["msp_name", "mod_id"]),
                          ("kegg_as_genes_sum", ["annotation"])]:
        files = merging_profiles.find_files_to_merge(samples, f"{pattern}.tsv.xz")
        merged_df = merging_profiles.merge_df(files, keys)
        numeric = merged_df.drop(columns=keys).to_numpy()
        row_sums = np.nansum(numeric, axis=1)
        occurrence = np.nansum((numeric != 0) & ~np.isnan(numeric), axis=1)
        expected = merged_df.loc[(row_sums >= merging_profiles.min_msp_abundance)
                                 & (occurrence >= merging_profiles.min_msp_occurrence), :]
        if remove_samples:
            expected = expected.loc[:, (expected.sum(axis=0) != 0)]
        expected.to_csv(tmp_path / f"{pattern}_expected.tsv", sep="\t", index=False)
        kept = merging_profiles._merge_filter_write(files, keys, tmp_path / f"{pattern}_fast.tsv")
        assert kept is not None
        assert (tmp_path / f"{pattern}_fast.tsv").read_text() == (
            tmp_path / f"{pattern}_expected.tsv"
        ).read_text(), pattern
        assert kept.equals(expected[keys].reset_index(drop=True)), pattern


def _write_profile(path: Path, rows: list[tuple]) -> Path:
    import lzma
    with lzma.open(path, "wt") as out:
        out.write("annotation\tvalue\n")
        for key, value in rows:
            out.write(f"{key}\t{value}\n")
    return path


def test_merge_filter_write_edge_cases(merging_profiles: Merging, tmp_path: Path, monkeypatch) -> None:
    """Cases handed back to the DataFrame path (None) and the pd.concat fallback"""
    merging_profiles.min_msp_occurrence = 1
    a = _write_profile(tmp_path / "a.tsv.xz", [("K1", 1.0), ("K2", 0.0)])
    b = _write_profile(tmp_path / "b.tsv.xz", [("K3", 2.5), ("K1", 0.5)])
    out = tmp_path / "out.tsv"
    keys = merging_profiles._merge_filter_write({"s1": a, "s2": b}, ["annotation"], out)
    assert list(keys["annotation"]) == ["K1", "K3"]
    expected = out.read_text()
    assert expected == "annotation\ts1\ts2\nK1\t1.0\t0.5\nK3\t\t2.5\n"
    # pandas helper failing (e.g. another signature): indexes combined with
    # pd.concat, same table. Only the first call (meteor's) fails.
    import pandas.core.indexes.api as pandas_api
    original = pandas_api._get_combined_index
    calls = []

    def first_call_fails(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise TypeError("unexpected signature")
        return original(*args, **kwargs)

    monkeypatch.setattr(pandas_api, "_get_combined_index", first_call_fails)
    assert merging_profiles._merge_filter_write({"s1": a, "s2": b}, ["annotation"], out) is not None
    assert len(calls) > 1
    assert out.read_text() == expected
    monkeypatch.undo()
    # duplicated key in a sample, missing key, quote in a sample name, nothing kept
    dup = _write_profile(tmp_path / "dup.tsv.xz", [("K1", 1.0), ("K1", 2.0)])
    assert merging_profiles._merge_filter_write({"s1": a, "s2": dup}, ["annotation"], out) is None
    missing = _write_profile(tmp_path / "nan.tsv.xz", [("", 1.0)])
    assert merging_profiles._merge_filter_write({"s1": missing}, ["annotation"], out) is None
    assert merging_profiles._merge_filter_write({'s"1': a}, ["annotation"], out) is None
    merging_profiles.min_msp_abundance = 1e9
    assert merging_profiles._merge_filter_write({"s1": a, "s2": b}, ["annotation"], out) is None


def test_execute_nothing_kept(merging_fast: Merging, tmp_path: Path) -> None:
    """No row passes the filters: tables written by the DataFrame path"""
    merging_fast.min_msp_abundance = 1e18
    merging_fast.execute()
    table = pd.read_table(tmp_path / "my_test_kegg_as_genes_sum.tsv")
    assert table.empty
    assert (tmp_path / "my_test_modules_completeness.tsv").exists()
