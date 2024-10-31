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

"""Test counter builder main objects"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..counter import Counter
from pathlib import Path
from hashlib import md5
from pysam import AlignmentFile
from itertools import chain
import pytest
import pandas as pd

# No more best count
# @pytest.fixture
# def counter_best(datadir: Path, tmp_path: Path) -> Counter:
#     meteor = Component
#     meteor.tmp_path = tmp_path
#     meteor.fastq_dir = datadir / "part1"
#     meteor.ref_dir = datadir / "catalogue"
#     meteor.ref_name = "test"
#     meteor.threads = 1
#     meteor.mapping_dir = tmp_path
#     return Counter(meteor, counting_type="best", mapping_type="end-to-end", trim=80,
#                    identity_threshold=0.95, alignment_number=1, counting_only=False,
#                    mapping_only=False)


@pytest.fixture
def counter_unique(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue" / "mock"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(
        meteor,
        counting_type="unique",
        mapping_type="end-to-end",
        trim=80,
        identity_user=0.95,
        alignment_number=10000,
        keep_filtered_alignments=True,
        identity_threshold=0.95,
        core_size=100,
    )


@pytest.fixture
def counter_smart_shared(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue" / "mock"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(
        meteor,
        counting_type="smart_shared",
        mapping_type="end-to-end",
        trim=80,
        identity_user=0.95,
        alignment_number=10000,
        keep_filtered_alignments=True,
        identity_threshold=0.95,
        core_size=100,
    )


@pytest.fixture
def counter_total(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue" / "mock"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(
        meteor,
        counting_type="total",
        mapping_type="end-to-end",
        trim=80,
        identity_user=0.95,
        alignment_number=10000,
        keep_all_alignments=False,
        keep_filtered_alignments=True,
        identity_threshold=0.95,
        core_size=100,
    )


def test_launch_mapping(counter_total: Counter):
    ref_json = counter_total.read_json(
        counter_total.meteor.ref_dir / "mock_reference.json"
    )
    census_json_file = counter_total.meteor.fastq_dir / "part1_census_stage_0.json"
    census_json = counter_total.read_json(census_json_file)
    sample_info = census_json["sample_info"]
    stage1_dir = counter_total.meteor.mapping_dir / sample_info["sample_name"]
    stage1_dir.mkdir(exist_ok=True, parents=True)
    counter_total.json_data[census_json_file] = {
        "census": census_json,
        "directory": stage1_dir,
        "Stage1FileName": stage1_dir
        / census_json_file.name.replace("stage_0", "stage_1"),
        "reference": ref_json,
    }
    counter_total.launch_mapping()
    assert counter_total.json_data[census_json_file]["Stage1FileName"].exists()
    # Fail with changing day
    # with counter_best.json_data[census_ini_file]["Stage1FileName"].open("rb") as stage1:
    #    assert md5(stage1.read()).hexdigest() == "a8a5b5e400dafb226ce3bab1a2cee69d"
    cram = stage1_dir / "part1_raw.cram"
    assert cram.exists()
    # cram = stage1_dir / "part1.cram"
    # bai = stage1_dir / "part1.cram.bai"
    # assert cram.exists()
    # assert bai.exists()


def test_write_table(counter_total: Counter, datadir: Path, tmp_path: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    output = tmp_path / "total.tsv.xz"
    counter_total.write_table(cramfile, output)
    assert output.exists()
    with output.open("rb") as out:
        assert md5(out.read()).hexdigest() == "006c2be62428ecf58cc20310d36078c1"


# detail of the mapping
# 2000 reads; of these:
#   2000 (100.00%) were unpaired; of these:
#     497 (24.85%) aligned 0 times
#     1454 (72.70%) aligned exactly 1 time
#     49 (2.45%) aligned >1 times
# 75.15% overall alignment rate
def test_filter_alignments(counter_total: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_total.filter_alignments(cramdesc)
        # We check that genes and reads are correctly associated
        # 11003 is mapped by two reads
        assert 26485 in genes["1368"]
        genes_list = list(set(map(int, chain.from_iterable(genes.values()))))
        # Detect correct reads within
        # read with 91.35802469135802 identity
        assert "1787" not in reads
        # read with multiple alignment with same identity 1.0 against gene 26485, 24457, 23617 and 24600
        assert "1368" in reads
        # 465 genes highlighted
        # All alignments are kept
        assert len(genes_list) == 9853


def test_uniq_from_mult(counter_unique: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_unique.filter_alignments(cramdesc)
        # print(genes)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult, unique_on_gene) = counter_unique.uniq_from_mult(
            reads, genes, database
        )
        genes_list = list(set(map(int, chain.from_iterable(genes_mult.values()))))
        assert len(genes_list) == 126
        assert "382" in unique_reads
        assert "1791" not in unique_reads
        assert "382" not in genes_mult
        assert "1791" in genes_mult
        assert 16622 in unique_on_gene


def test_compute_co(counter_smart_shared: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_smart_shared.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (
            _,
            genes_mult,
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, co = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        # We check genes with no unique counts, but multiple reads
        assert co[("1368", 26485)] == 0.25
        assert "1368" in read_dict[26485]
        # We check multiple alignment on the same read
        assert co[("11619", 6205)] == 1.0
        # We check the normal case
        assert co[("11450", 18783)] == 2 / 3


def test_get_co_coefficient(counter_smart_shared: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_smart_shared.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (
            _,
            genes_mult,  # pylint: disable=unused-variable
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, co_dict = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        assert (
            next(counter_smart_shared.get_co_coefficient(18783, read_dict, co_dict))
            == 2 / 3
        )


def test_compute_abm(counter_smart_shared: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_smart_shared.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (
            _,
            genes_mult,  # pylint: disable=unused-variable
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, coef_read = counter_smart_shared.compute_co(
            genes_mult, unique_on_gene
        )
        multiple_dict = counter_smart_shared.compute_abm(read_dict, coef_read, database)
        assert multiple_dict[18783] == 2 / 3


def test_compute_abs(counter_smart_shared: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_smart_shared.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (
            _,
            genes_mult,  # pylint: disable=unused-variable
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, coef_read = counter_smart_shared.compute_co(
            genes_mult, unique_on_gene
        )
        multiple_dict = counter_smart_shared.compute_abm(read_dict, coef_read, database)
        abundance = counter_smart_shared.compute_abs(unique_on_gene, multiple_dict)
        assert abundance[18783] == 2 + 2 / 3


def test_write_stat(counter_smart_shared: Counter, tmp_path: Path) -> None:
    test = tmp_path / "test.tsv.xz"
    abundance_dict = {"11003": 1.5}
    database = {"11003": 500}
    counter_smart_shared.write_stat(test, abundance_dict, database)
    with test.open("rb") as out:
        assert md5(out.read()).hexdigest() == "64cc59536603837848cbbfbb563d04fc"


def test_save_cram(counter_unique: Counter, datadir: Path, tmp_path: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    ref_json = counter_unique.read_json(
        counter_unique.meteor.ref_dir / "mock_reference.json"
    )
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, _ = counter_unique.filter_alignments(
            cramdesc
        )  # pylint: disable=unused-variable
        read_list = reads.values()
        merged_list = chain.from_iterable(read_list)
        tmpcramfile = tmp_path / "test"
        counter_unique.save_cram_strain(tmpcramfile, cramdesc, merged_list, ref_json)
        assert tmpcramfile.exists()
        # issues at testing content
        # with tmpcramfile.open("rb") as out:
        #     assert md5(out.read()).hexdigest() == "ef47276b6fcb7ad398801b7f5c52ef04"


def test_launch_counting_unique(counter_unique: Counter, datadir: Path, tmp_path: Path):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = counter_unique.meteor.fastq_dir / "part1_census_stage_0.json"
    census_json = counter_unique.read_json(census_json_file)
    ref_json = counter_unique.read_json(
        counter_unique.meteor.ref_dir / "mock_reference.json"
    )
    counter_unique.launch_counting(
        raw_cramfile, cramfile, countfile, ref_json, census_json, census_json_file
    )
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "f5bc528dcbf594b5089ad7f6228ebab5"


def test_launch_counting_total(counter_total: Counter, datadir: Path, tmp_path: Path):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = counter_total.meteor.fastq_dir / "part1_census_stage_0.json"
    census_json = counter_total.read_json(census_json_file)
    ref_json = counter_total.read_json(
        counter_total.meteor.ref_dir / "mock_reference.json"
    )
    counter_total.launch_counting(
        raw_cramfile, cramfile, countfile, ref_json, census_json, census_json_file
    )
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "f010e4136323ac408d4c127e243756c2"


def test_launch_counting_smart_shared(
    counter_smart_shared: Counter, datadir: Path, tmp_path: Path
):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = (
        counter_smart_shared.meteor.fastq_dir / "part1_census_stage_0.json"
    )
    census_json = counter_smart_shared.read_json(census_json_file)
    ref_json = counter_smart_shared.read_json(
        counter_smart_shared.meteor.ref_dir / "mock_reference.json"
    )
    counter_smart_shared.launch_counting(
        raw_cramfile, cramfile, countfile, ref_json, census_json, census_json_file
    )
    # with countfile.open("rb") as out:
    #     assert md5(out.read()).hexdigest() == "4bdd7327cbad8e71d210feb0c6375077"
    expected_output = pd.read_csv(
        datadir / "expected_output" / "count.tsv.xz", sep="\t"
    )
    count_data = pd.read_csv(countfile, sep="\t")
    assert count_data.equals(expected_output)


def test_execute(counter_smart_shared: Counter, tmp_path: Path):
    counter_smart_shared.execute()
    part1 = tmp_path / "part1/part1.tsv.xz"
    assert part1.exists()
    with part1.open("rb") as out:
        assert md5(out.read()).hexdigest() == "5db950a4404793f73ba034e99cb676fa"
