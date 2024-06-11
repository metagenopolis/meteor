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

"""Test profiler main object"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..profiler import Profiler
from pathlib import Path
import pytest
import pandas as pd


@pytest.fixture
def profiler_standard(datadir: Path, tmp_path: Path) -> Profiler:
    meteor = Component
    meteor.mapping_dir = datadir / "mapping"
    meteor.profile_dir = tmp_path
    meteor.ref_dir = datadir / "catalogue"
    return Profiler(
        meteor=meteor,
        rarefaction_level=-1,
        seed=12345,
        coverage_factor=100.0,
        normalization=None,
        core_size=4,
        msp_filter_user=0.5,
        completeness=0.6,
        msp_filter=0.5,
    )


def test_rarefy(profiler_standard: Profiler, datadir: Path) -> None:
    # No rarefaction with 0 unmapped reads
    profiler_standard.rarefy(rarefaction_level=900, unmapped_reads=0, seed=12345)
    expected_output = pd.read_table(datadir / "expected_output" / "sample_rarefied.tsv")
    assert profiler_standard.gene_count.equals(expected_output)


def test_rarefy2(profiler_standard: Profiler, datadir: Path) -> None:
    # Rarefaction with 0 unmapped reads
    profiler_standard.rarefy(rarefaction_level=100, unmapped_reads=0, seed=12345)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_raref_100_unmapped_0.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_rarefy3(profiler_standard: Profiler, datadir: Path) -> None:
    # Rarefaction with 100 unmapped reads
    profiler_standard.rarefy(rarefaction_level=900, unmapped_reads=100, seed=12345)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_raref_900_unmapped_100.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_rarefy4(profiler_standard: Profiler, datadir: Path) -> None:
    # No rarefaction with unmapped reads
    profiler_standard.rarefy(rarefaction_level=900, unmapped_reads=80, seed=12345)
    expected_output = pd.read_table(datadir / "expected_output" / "sample_rarefied.tsv")
    assert profiler_standard.gene_count.equals(expected_output)


def test_normalize_coverage(profiler_standard: Profiler, datadir: Path) -> None:
    # Manually change gene count so that all possibilities are tested
    manual_gene_count = {
        "gene_id": [1, 2, 3, 4, 5, 6],
        "gene_length": [80, 81, 160, 159, 200, 2000],
        "value": [27, 58, 12, 99, 80, 89],
    }
    profiler_standard.gene_count = pd.DataFrame(manual_gene_count)
    profiler_standard.normalize_coverage(trim_length=80)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_coverage.tsv"
    )
    assert profiler_standard.gene_count.round(6).equals(expected_output.round(6))


def test_normalize_fpkm1(profiler_standard: Profiler, datadir: Path) -> None:
    # No downsizing, no unmapped reads
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=0)
    # Round the results to 6 digits, if not df are not equal
    profiler_standard.gene_count["value"] = (
        profiler_standard.gene_count["value"].astype(float).round(6)
    )

    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_fpkm_raref_0_unmapped_0.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_normalize_fpkm2(profiler_standard: Profiler, datadir: Path) -> None:
    # No downsizing, 500 unmapped reads
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=500)
    # Round the results to 6 digits, if not df are not equal
    profiler_standard.gene_count["value"] = (
        profiler_standard.gene_count["value"].astype(float).round(6)
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_fpkm_raref_0_unmapped_500.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_normalize_fpkm3(profiler_standard: Profiler, datadir: Path) -> None:
    # Downsizing at 100, no unmapped reads
    profiler_standard.rarefy(rarefaction_level=100, unmapped_reads=0, seed=12345)
    profiler_standard.normalize_fpkm(rarefaction_level=100, unmapped_reads=0)
    # Round the results to 6 digits, if not df are not equal
    profiler_standard.gene_count["value"] = profiler_standard.gene_count["value"].round(
        6
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_fpkm_raref_100_unmapped_0.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_normalize_fpkm4(profiler_standard: Profiler, datadir: Path) -> None:
    # Downsizing at 900, 100 unmapped reads
    profiler_standard.rarefy(rarefaction_level=900, unmapped_reads=100, seed=12345)
    profiler_standard.normalize_fpkm(rarefaction_level=900, unmapped_reads=100)
    # Round the results to 6 digits, if not df are not equal
    profiler_standard.gene_count["value"] = profiler_standard.gene_count["value"].round(
        6
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_fpkm_raref_900_unmapped_100.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_normalize_fpkm5(profiler_standard: Profiler, datadir: Path) -> None:
    # Downsizing at 900, 80 unmapped reads => no downsizing (890 total reads)
    profiler_standard.rarefy(rarefaction_level=900, unmapped_reads=80, seed=12345)
    profiler_standard.normalize_fpkm(rarefaction_level=900, unmapped_reads=80)
    # Round the results to 6 digits, if not df are not equal
    profiler_standard.gene_count["value"] = profiler_standard.gene_count["value"].round(
        6
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_norm_fpkm_raref_900_unmapped_80.tsv"
    )
    assert profiler_standard.gene_count.equals(expected_output)


def test_compute_msp(profiler_standard: Profiler, datadir: Path) -> None:
    # Define msp core genes
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    # Compute MSP directly on the raw gene count (only rounded)
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    # Round at 6 digits
    profiler_standard.msp_table["value"] = (
        profiler_standard.msp_table["value"].astype(float).round(6)
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_msp.tsv"
    )
    assert profiler_standard.msp_table.equals(expected_output)
    # Compute MSP on fpkm gene count
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=0)
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    # Round at 6 digits
    profiler_standard.msp_table["value"] = (
        profiler_standard.msp_table["value"].astype(float).round(6)
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_norm_fpkm_msp.tsv"
    )
    assert profiler_standard.msp_table.equals(expected_output)


def test_get_msp_core(profiler_standard: Profiler) -> None:
    # Define msp core genes
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    true_set = {
        "msp_01_standard": {5, 11, 17, 18},
        "msp_02_undetected": {1, 10, 13, 15},
        "msp_03_small": {10, 13, 14},
        "msp_04_small_nd": {3, 8, 19},
    }
    assert msp_set == true_set


def test_compute_msp_stats(profiler_standard: Profiler) -> None:
    # Compute stats on the raw gene matrix
    msp_stats_raw_matrix = profiler_standard.compute_msp_stats(
        msp_def_filename=profiler_standard.msp_filename
    )
    assert msp_stats_raw_matrix == 0.94
    # Compute stats on the fpkm gene matrix
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=0)
    msp_stats_matrix_fpkm = profiler_standard.compute_msp_stats(
        msp_def_filename=profiler_standard.msp_filename
    )
    assert msp_stats_matrix_fpkm == 0.97


def test_compute_ko_abundance(profiler_standard: Profiler, datadir: Path) -> None:
    # Compute KO on the raw gene matrix
    profiler_standard.compute_ko_abundance(
        annot_file=profiler_standard.db_filenames["kegg"]
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_kegg.tsv"
    )
    profiler_standard.functions["value"] = (
        profiler_standard.functions["value"].astype(float).round(6)
    )
    assert profiler_standard.functions.equals(expected_output)
    # Compute KO on the fpkm gene matrix
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=0)
    profiler_standard.compute_ko_abundance(
        annot_file=profiler_standard.db_filenames["kegg"]
    )
    # Round at 6 digits
    profiler_standard.functions["value"] = (
        profiler_standard.functions["value"].astype(float).round(6)
    )
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_norm_fpkm_kegg.tsv"
    )
    assert profiler_standard.functions.equals(expected_output)


def test_compute_ko_abundance_by_msp(
    profiler_standard: Profiler, datadir: Path
) -> None:
    # Compute KO on the raw gene matrix
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    profiler_standard.compute_ko_abundance_by_msp(
        annot_file=profiler_standard.db_filenames["kegg"],
        msp_def_filename=profiler_standard.msp_filename,
    )
    # Round at 6 digits
    profiler_standard.functions["value"] = profiler_standard.functions["value"].round(6)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_kegg_by_msp.tsv"
    )
    assert profiler_standard.functions.equals(expected_output)
    # Compute KO on the fpkm gene matrix
    profiler_standard.normalize_fpkm(rarefaction_level=0, unmapped_reads=0)
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    profiler_standard.compute_ko_abundance_by_msp(
        annot_file=profiler_standard.db_filenames["kegg"],
        msp_def_filename=profiler_standard.msp_filename,
    )
    # Round at 6 digits
    profiler_standard.functions["value"] = profiler_standard.functions["value"].round(6)

    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_norm_fpkm_kegg_by_msp.tsv"
    )
    assert profiler_standard.functions.equals(expected_output)


def test_compute_ko_stats(profiler_standard: Profiler) -> None:
    # Compute stats on ko computed as sum of genes
    # Use the raw gene cout only
    ko_stats = profiler_standard.compute_ko_stats(
        annot_file=profiler_standard.db_filenames["kegg"],
        by_msp=False,
        msp_def_filename=profiler_standard.msp_filename,
    )
    assert ko_stats == 0.83
    # Compute stats on ko computed as sum of MSP
    # Use the raw gene cout only
    ko_stats = profiler_standard.compute_ko_stats(
        annot_file=profiler_standard.db_filenames["kegg"],
        by_msp=True,
        msp_def_filename=profiler_standard.msp_filename,
    )
    assert ko_stats == 0.82


def test_merge_catalogue_info(profiler_standard: Profiler, datadir: Path) -> None:
    # Compute MSP
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    # Merge information for KEGG and MUSTARD
    db_to_merge = {db: profiler_standard.db_filenames[db] for db in ["kegg", "mustard"]}
    merged_catalogue = profiler_standard.merge_catalogue_info(
        msp_file=profiler_standard.msp_filename,
        annot_file=db_to_merge,
    )
    # Order rows
    merged_catalogue = merged_catalogue.sort_values(
        by=[
            "msp_name",
            "gene_id",
            "annotation",
        ]
    ).reset_index(drop=True)
    expected_output = pd.read_table(
        datadir / "expected_output" / "merged_catalogue_info.tsv"
    )
    expected_output = expected_output.sort_values(
        by=[
            "msp_name",
            "gene_id",
            "annotation",
        ]
    ).reset_index(drop=True)
    assert merged_catalogue.equals(expected_output)


def test_compute_completeness(profiler_standard: Profiler) -> None:
    # Compute MGS table
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    # Annotate gene from MSP with KO
    annotated_msp = profiler_standard.merge_catalogue_info(
        profiler_standard.msp_filename, profiler_standard.db_filenames
    )
    module_alternatives = [{"K8", "K9"}, {"TIGRFAM01", "K9"}, {"K8", "ARD1", "ARD2"}]
    completeness = profiler_standard.compute_completeness(
        mod=module_alternatives, annotated_msp=annotated_msp
    )
    # Round the result to 2 digits
    completeness = {key: round(value, 2) for key, value in completeness.items()}
    expected_output = {
        "msp_01_standard": 0.67,
        "msp_02_undetected": 0.5,
        "msp_03_small": 0.5,
        "msp_04_small_nd": 0.33,
    }
    assert completeness == expected_output


def test_compute_max(profiler_standard: Profiler) -> None:
    # With a complete alternative
    module_alternatives = [{"ARD1"}, {"K02", "KO1", "ARD1", "K03"}, {"ARD3", "ARD4"}]
    set_ko = {"K01", "K02", "ARD1", "ARD2"}
    max_completeness = profiler_standard.compute_max(mod=module_alternatives, ko=set_ko)
    assert max_completeness == 1.0
    # With zero common KO
    module_alternatives = [{"ARD3", "ARD4"}]
    set_ko = {"K01", "K02", "ARD1", "ARD2"}
    max_completeness = profiler_standard.compute_max(mod=module_alternatives, ko=set_ko)
    assert max_completeness == 0.0
    # With a fraction
    module_alternatives = [
        {"ARD3", "ARD4"},
        {"K02", "K03", "K04"},
        {"K01", "K03", "K04"},
    ]
    set_ko = {"K01", "K02", "ARD1", "ARD2"}
    max_completeness = profiler_standard.compute_max(mod=module_alternatives, ko=set_ko)
    assert round(max_completeness, 2) == 0.33


def test_compute_completeness_all(profiler_standard: Profiler) -> None:
    # Compute MGS table
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    # Annotate gene from MSP with KO
    annotated_msp = profiler_standard.merge_catalogue_info(
        profiler_standard.msp_filename, profiler_standard.db_filenames
    )
    # Define list of modules
    all_mod = {
        "mod1": [{"ARD2"}, {"ARD1"}],
        "mod2": [{"K8"}, {"K9"}],
        "mod3": [{"K2", "K3", "K4"}, {"K1"}, {"K2", "K4"}, {"K1", "K8"}],
        "mod4": [{"fonction1", "fonction2"}],
    }
    # Compute completeness for all modules in all msp
    completeness = profiler_standard.compute_completeness_all(
        all_mod=all_mod, annotated_msp=annotated_msp
    )
    expected_output = {
        "mod1": {
            "msp_01_standard": 1.0,
            "msp_02_undetected": 0.0,
            "msp_03_small": 0.0,
            "msp_04_small_nd": 1.0,
        },
        "mod2": {
            "msp_01_standard": 1.0,
            "msp_02_undetected": 1.0,
            "msp_03_small": 1.0,
            "msp_04_small_nd": 0.0,
        },
        "mod3": {
            "msp_01_standard": 1.0,
            "msp_02_undetected": 0.5,
            "msp_03_small": 0.5,
            "msp_04_small_nd": 0.5,
        },
        "mod4": {
            "msp_01_standard": 0.0,
            "msp_02_undetected": 0.0,
            "msp_03_small": 0.0,
            "msp_04_small_nd": 0.0,
        },
    }
    assert completeness == expected_output


def test_compute_module_abundance(profiler_standard: Profiler, datadir: Path) -> None:
    # Compute MGS table
    msp_set = profiler_standard.get_msp_core(
        profiler_standard.msp_filename, core_size=4
    )
    profiler_standard.compute_msp(msp_dict=msp_set, filter_pc=0.5)
    all_mod = {
        "mod1": [{"ARD2"}, {"ARD1"}],
        "mod2": [{"K8"}, {"K9"}],
        "mod3": [{"K2", "K3", "K4"}, {"K1"}, {"K2", "K4"}, {"K1", "K8"}],
        "mod4": [{"fonction1", "fonction2"}],
    }
    profiler_standard.compute_module_abundance(
        msp_file=profiler_standard.msp_filename,
        annot_file=profiler_standard.db_filenames,
        all_mod=all_mod,
        completeness=0.6,
    )
    # Round at 6 digits
    profiler_standard.mod_table["value"] = profiler_standard.mod_table["value"].round(6)
    profiler_standard.mod_completeness["value"] = profiler_standard.mod_completeness[
        "value"
    ].round(6)
    # Check module abundance matrix
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_modules.tsv"
    )
    assert profiler_standard.mod_table.equals(expected_output)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_mod_completeness.tsv"
    )
    # Check module completeness matrix
    assert profiler_standard.mod_completeness.equals(expected_output)


def test_execute(profiler_standard: Profiler, datadir: Path) -> None:
    profiler_standard.module_path = datadir / "module.tsv"
    profiler_standard.execute()
    # Check symlink file (raw data)
    raw_gene_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_raw.tsv.xz"
    )
    assert raw_gene_table_file.exists()
    assert raw_gene_table_file.is_symlink()
    expected_output = pd.read_table(datadir / "mapping" / "sample.tsv.xz")
    assert pd.read_table(raw_gene_table_file).equals(expected_output)
    # Check gene file
    gene_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_genes.tsv.xz"
    )
    gene_table = pd.read_table(gene_table_file, compression="xz")
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rarefaction.tsv"
    )
    assert gene_table.equals(expected_output)
    # Check MSP file
    msp_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_msp.tsv.xz"
    )
    msp_table = pd.read_table(msp_table_file, compression="xz")
    msp_table["value"] = msp_table["value"].round(6)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_msp.tsv"
    )
    assert msp_table.equals(expected_output)
    # Check functions as sum of genes (mustard only)
    fun_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_mustard_as_genes_sum.tsv.xz"
    )
    fun_table = pd.read_table(fun_table_file, compression="xz")
    fun_table["value"] = fun_table["value"].round(6)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_mustard.tsv"
    )
    assert fun_table.equals(expected_output)
    # Check functions as sum of MSP
    fun_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_dbcan_as_msp_sum.tsv.xz"
    )
    fun_table = pd.read_table(fun_table_file, compression="xz")
    fun_table["value"] = fun_table["value"].round(6)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_dbcan_by_msp.tsv"
    )
    assert fun_table.equals(expected_output)
    # Check modules
    # 1. Module abundance table
    module_table_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_modules.tsv.xz"
    )
    module_table = pd.read_table(module_table_file, compression="xz")
    module_table["value"] = module_table["value"].round(6)
    expected_output = pd.read_table(
        datadir / "expected_output" / "sample_no_rf_no_norm_modules_kegg_only.tsv"
    )
    assert module_table.equals(expected_output)
    # 2. Module completeness table
    module_completeness_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_modules_completeness.tsv.xz"
    )
    module_completeness = pd.read_table(module_completeness_file, compression="xz")
    module_completeness["value"] = module_completeness["value"].round(6)
    expected_output = pd.read_table(
        datadir
        / "expected_output"
        / "sample_no_rf_no_norm_mod_completeness_kegg_only.tsv"
    )
    assert module_completeness.equals(expected_output)
    census_stage_2_file = (
        profiler_standard.stage2_dir
        / f"{profiler_standard.output_base_filename}_census_stage_2.json"
    )
    assert census_stage_2_file.exists()
