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
from configparser import ConfigParser
from pysam import AlignmentFile
from itertools import chain
import pytest

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
        identity_threshold=0.95,
        alignment_number=10000,
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
        identity_threshold=0.95,
        alignment_number=10000,
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
        identity_threshold=0.95,
        alignment_number=10000,
    )


def test_launch_mapping(counter_total: Counter):
    ref_ini_file = counter_total.meteor.ref_dir / "mock_reference.ini"
    ref_ini = ConfigParser()
    with open(ref_ini_file, "rt", encoding="UTF-8") as ref:
        ref_ini.read_file(ref)
    census_ini_file = counter_total.meteor.fastq_dir / "part1_census_stage_0.ini"
    census_ini = ConfigParser()
    with open(census_ini_file, "rt", encoding="UTF-8") as cens:
        census_ini.read_file(cens)
    sample_info = census_ini["sample_info"]
    stage1_dir = counter_total.meteor.mapping_dir / sample_info["sample_name"]
    stage1_dir.mkdir(exist_ok=True, parents=True)
    counter_total.ini_data[census_ini_file] = {
        "census": census_ini,
        "directory": stage1_dir,
        "Stage1FileName": stage1_dir
        / census_ini_file.name.replace("stage_0", "stage_1"),
        "reference": ref_ini,
    }
    counter_total.launch_mapping()
    assert counter_total.ini_data[census_ini_file]["Stage1FileName"].exists()
    # Fail with changing day
    # with counter_best.ini_data[census_ini_file]["Stage1FileName"].open("rb") as stage1:
    #    assert md5(stage1.read()).hexdigest() == "a8a5b5e400dafb226ce3bab1a2cee69d"
    sam = stage1_dir / "part1.sam"
    assert sam.exists()
    # bam = stage1_dir / "part1.bam"
    # bai = stage1_dir / "part1.bam.bai"
    # assert bam.exists()
    # assert bai.exists()


def test_write_table(counter_total: Counter, datadir: Path, tmp_path: Path) -> None:
    bamfile = datadir / "total.bam"
    output = tmp_path / "total.tsv"
    counter_total.write_table(bamfile, output)
    assert output.exists()
    with output.open("rb") as out:
        assert md5(out.read()).hexdigest() == "647decb83957c2d3257c3a13fd399c92"


# detail of the mapping
# 2000 reads; of these:
#   2000 (100.00%) were unpaired; of these:
#     497 (24.85%) aligned 0 times
#     1454 (72.70%) aligned exactly 1 time
#     49 (2.45%) aligned >1 times
# 75.15% overall alignment rate
def test_filter_alignments(counter_total: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_total.filter_alignments(bamdesc)
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
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_unique.filter_alignments(bamdesc)
        # print(genes)
        references = map(int, bamdesc.references)
        # get reference length
        lengths = bamdesc.lengths
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
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_alignments(bamdesc)
        references = map(int, bamdesc.references)
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (
            unique_reads,
            genes_mult,  # pylint: disable=unused-variable
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
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_alignments(bamdesc)
        references = map(int, bamdesc.references)
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (
            unique_reads,
            genes_mult,  # pylint: disable=unused-variable
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, co_dict = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        assert (
            next(counter_smart_shared.get_co_coefficient(18783, read_dict, co_dict))
            == 2 / 3
        )


def test_compute_abm(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_alignments(bamdesc)
        references = map(int, bamdesc.references)
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (
            unique_reads,
            genes_mult,  # pylint: disable=unused-variable
            unique_on_gene,
        ) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, coef_read = counter_smart_shared.compute_co(
            genes_mult, unique_on_gene
        )
        multiple_dict = counter_smart_shared.compute_abm(read_dict, coef_read, database)
        assert multiple_dict[18783] == 2 / 3


def test_compute_abs(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_alignments(bamdesc)
        references = map(int, bamdesc.references)
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (
            unique_reads,
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
    test = tmp_path / "test.tsv"
    abundance_dict = {"11003": 1.5}
    database = {"11003": 500}
    counter_smart_shared.write_stat(test, abundance_dict, database)
    with test.open("rb") as out:
        assert md5(out.read()).hexdigest() == "0448f393b702b038840f1be20c0f4aa6"


def test_save_bam(counter_unique: Counter, datadir: Path, tmp_path: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_unique.filter_alignments(
            bamdesc
        )  # pylint: disable=unused-variable
        read_list = list(chain(reads.values()))
        merged_list = list(chain.from_iterable(read_list))
        tmpbamfile = tmp_path / "test"
        counter_unique.save_bam(tmpbamfile, bamdesc, merged_list)
        assert tmpbamfile.exists()
        # issues at testing content
        # with tmpbamfile.open("rb") as out:
        #     assert md5(out.read()).hexdigest() == "7e9c1b3e89690624ca03882cb968fb09"


def test_launch_counting_unique(counter_unique: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_unique.launch_counting(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "4188e42b16ca23af21bb0b03c88089fe"


def test_launch_counting_total(counter_total: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_total.launch_counting(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "78216b228518d101983f39a92ac9bdb0"


def test_launch_counting_smart_shared(
    counter_smart_shared: Counter, datadir: Path, tmp_path: Path
):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_smart_shared.launch_counting(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "06383f066a83fe158f9d6eb98c9a41a3"


def test_execute(counter_smart_shared: Counter, tmp_path: Path):
    counter_smart_shared.execute()
    part1 = tmp_path / "part1/part1.tsv"
    assert part1.exists()
    with part1.open("rb") as out:
        assert md5(out.read()).hexdigest() == "91ee131f2cbbb3abdb441849b55c4481"
