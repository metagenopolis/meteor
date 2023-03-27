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


@pytest.fixture
def counter_best(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(meteor, counting_type="best", mapping_type="end-to-end", trim=80,
                   identity_threshold=0.95, alignment_number=1, counting_only=False,
                   mapping_only=False)


@pytest.fixture
def counter_unique(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(meteor, counting_type="unique", mapping_type="end-to-end", trim=80,
                   identity_threshold=0.95, alignment_number=10000, counting_only=False, mapping_only=False)


@pytest.fixture
def counter_smart_shared(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(meteor, counting_type="smart_shared", mapping_type="end-to-end", trim=80,
                   identity_threshold=0.95, alignment_number=10000, counting_only=False, mapping_only=False)


@pytest.fixture
def counter_total(datadir: Path, tmp_path: Path) -> Counter:
    meteor = Component
    meteor.tmp_path = tmp_path
    meteor.fastq_dir = datadir / "part1"
    meteor.ref_dir = datadir / "catalogue"
    meteor.ref_name = "test"
    meteor.threads = 1
    meteor.mapping_dir = tmp_path
    return Counter(meteor, counting_type="total", mapping_type="end-to-end", trim=80,
                   identity_threshold=0.95, alignment_number=10000, counting_only=False, mapping_only=False)


def test_launch_mapping2(counter_best: Counter):
    ref_ini_file = counter_best.meteor.ref_dir / "mock" / "mock_reference.ini"
    ref_ini = ConfigParser()
    with open(ref_ini_file, "rt", encoding="UTF-8") as ref:
        ref_ini.read_file(ref)
    census_ini_file = counter_best.meteor.fastq_dir / "part1_census_stage_0.ini"
    census_ini = ConfigParser()
    with open(census_ini_file, "rt", encoding="UTF-8") as cens:
        census_ini.read_file(cens)
    sample_info = census_ini["sample_info"]
    stage1_dir = counter_best.meteor.mapping_dir / sample_info["sample_name"]
    stage1_dir.mkdir(exist_ok=True, parents=True)
    counter_best.ini_data[census_ini_file] = {
        "census": census_ini,
        "directory": stage1_dir,
        "Stage1FileName": stage1_dir / census_ini_file.name.replace("stage_0", "stage_1"),
        "reference": ref_ini
    }
    counter_best.launch_mapping2()
    assert counter_best.ini_data[census_ini_file]["Stage1FileName"].exists()
    # Fail with changing day
    # with counter_best.ini_data[census_ini_file]["Stage1FileName"].open("rb") as stage1:
    #    assert md5(stage1.read()).hexdigest() == "a8a5b5e400dafb226ce3bab1a2cee69d"
    bam = stage1_dir / "part1.bam"
    bai = stage1_dir / "part1.bam.bai"
    assert bam.exists()
    assert bai.exists()


def test_write_table(counter_best: Counter, datadir: Path, tmp_path: Path) -> None:
    bamfile = datadir / "best.bam"
    output = tmp_path / "best.tsv"
    counter_best.write_table(bamfile, output)
    assert output.exists()
    with output.open("rb") as out:
        assert md5(out.read()).hexdigest() == "b8b640ba131d9f84483f22571b50aa0b"


# detail of the mapping
# 500 reads; of these:
#   500 (100.00%) were unpaired; of these:
#     119 (23.80%) aligned 0 times
#     366 (73.20%) aligned exactly 1 time
#     15 (3.00%) aligned >1 times
def test_filter_bam_best(counter_best: Counter, datadir: Path) -> None:
    bamfile = datadir / "best.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_best.filter_bam(bamdesc)
        genes_list = list(set(map(int, chain.from_iterable(genes.values()))))
        # read with 93.75 identity
        assert "7ZVI4:08019:00760_2" not in reads
        # read with 96.25 identity
        assert "7ZVI4:08031:00789_2" in reads
        # 358 genes highlighted
        assert len(genes_list) == 358


def test_filter_bam_total(counter_total: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_total.filter_bam(bamdesc)
        # 11003 is mapped by two reads
        assert "11003" in genes["7ZVI4:02650:13780_2"]
        assert "11003" in genes["7ZVI4:02660:14094_2"]
        genes_list = list(set(map(int, chain.from_iterable(genes.values()))))
        # read with 93.75 identity
        assert "7ZVI4:02653:13846_2" not in reads
        # read with multiple alignment with same score (148) against gene 3971 and 5259
        assert "7ZVI4:02660:13822_2" in reads
        # read with multiple alignment with same score (0) against gene 11003 and 26832
        assert "7ZVI4:02650:13780_2" in reads
        # 465 genes highlighted
        # All alignments are kept
        assert len(genes_list) == 465


def test_filter_bam_smart_shared(counter_smart_shared: Counter, datadir: Path) -> None:
    # unique and shared counting have at this step the same result
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_bam(bamdesc)
        # 11003 is mapped by two reads
        assert "11003" in genes["7ZVI4:02650:13780_2"]
        assert "11003" in genes["7ZVI4:02660:14094_2"]
        genes_list = list(set(map(int, chain.from_iterable(genes.values()))))
        # read with 93.75 identity
        assert "7ZVI4:02653:13846_2" not in reads
        # read with multiple alignment with same score (148) against gene 3971 and 5259
        assert "7ZVI4:02660:13822_2" in reads
        # read with multiple alignment with same score (0) against gene 11003 and 26832
        assert "7ZVI4:02650:13780_2" in reads
        # 459 genes highlighted
        # All alignments are kept
        assert len(genes_list) == 459


def test_filter_bam_unique(counter_unique: Counter, datadir: Path) -> None:
    # Total and shared counting have at this step the same result
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_unique.filter_bam(bamdesc)
        genes_list = list(set(map(int, chain.from_iterable(genes.values()))))
        # read with 93.75 identity
        assert "7ZVI4:02653:13846_2" not in reads
        # read with multiple alignment with same score (148) against gene 3971 and 5259
        assert "7ZVI4:02660:13822_2" in reads
        # read with multiple alignment with same score (0) against gene 11003 and 26832
        assert "7ZVI4:02650:13780_2" in reads
        # 459 genes highlighted
        # All alignments are kept
        assert len(genes_list) == 459


def test_uniq_from_mult(counter_unique: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_unique.filter_bam(bamdesc)
        references = bamdesc.references
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads,
         genes_mult,
         unique_on_gene) = counter_unique.uniq_from_mult(reads, genes, database)  # pylint: disable=unused-variable
        assert "7ZVI4:02660:13822_2" not in unique_reads
        genes_list = list(set(map(int, chain.from_iterable(genes_mult.values()))))
        assert len(genes_list) == 59


def test_compute_co(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_bam(bamdesc)
        # 11003 is mapped by two reads
        assert "11003" in genes["7ZVI4:02650:13780_2"]
        assert "11003" in genes["7ZVI4:02660:14094_2"]
        references = bamdesc.references
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult,                             # pylint: disable=unused-variable
         unique_on_gene) = counter_smart_shared.uniq_from_mult(reads, genes,
                                                               database)
        read_dict, co = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        # gene "7ZVI4:02660:13822_2" did not align with gene aligned by other reads
        assert co[("7ZVI4:02650:13780_2", "11003")] == 0.5
        assert "7ZVI4:02650:13780_2" in read_dict["11003"]


def test_get_co_coefficient(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_bam(bamdesc)
        references = bamdesc.references
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult,                   # pylint: disable=unused-variable
         unique_on_gene) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, co_dict = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        assert next(counter_smart_shared.get_co_coefficient("11003", read_dict, co_dict)) == 0.5


def test_compute_abm(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_bam(bamdesc)
        references = bamdesc.references
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult,                 # pylint: disable=unused-variable
         unique_on_gene) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, coef_read = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        multiple_dict = counter_smart_shared.compute_abm(read_dict, coef_read, database)
        assert multiple_dict["11003"] == 0.5


def test_compute_abs(counter_smart_shared: Counter, datadir: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_smart_shared.filter_bam(bamdesc)
        references = bamdesc.references
        # get reference length
        lengths = bamdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult,                # pylint: disable=unused-variable
         unique_on_gene) = counter_smart_shared.uniq_from_mult(reads, genes, database)
        read_dict, coef_read = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
        multiple_dict = counter_smart_shared.compute_abm(read_dict, coef_read, database)
        abundance = counter_smart_shared.compute_abs(unique_on_gene, multiple_dict)
        assert abundance["11003"] == 1.5


def test_write_stat(counter_smart_shared: Counter, tmp_path: Path) -> None:
    test = tmp_path / "test.tsv"
    abundance_dict = {"11003": 1.5}
    database = {"11003": 500}
    counter_smart_shared.write_stat(test, abundance_dict, database)
    with test.open("rb") as out:
        assert md5(out.read()).hexdigest() == "0df6018a478d2c7d9d3b57897edbf60c"


def test_save_bam(counter_unique: Counter, datadir: Path, tmp_path: Path) -> None:
    bamfile = datadir / "total.bam"
    with AlignmentFile(str(bamfile.resolve()), "rb") as bamdesc:
        reads, genes = counter_unique.filter_bam(bamdesc)  # pylint: disable=unused-variable
        read_list = list(chain(reads.values()))
        merged_list = list(chain.from_iterable(read_list))
        tmpbamfile = tmp_path / "test"
        counter_unique.save_bam(tmpbamfile, bamdesc, merged_list)
        assert tmpbamfile.exists()
        # issues at testing content
        # with tmpbamfile.open("rb") as out:
        #     assert md5(out.read()).hexdigest() == "7e9c1b3e89690624ca03882cb968fb09"


def test_launch_counting2_unique(counter_unique: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_unique.launch_counting2(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "9e0d47640c24f5c31fa49eed171a04e3"


def test_launch_counting2_best(counter_best: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "best.bam"
    countfile = tmp_path / "count.tsv"
    counter_best.launch_counting2(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "8539ee3303d283bd74af3d11a6c284bb"


def test_launch_counting2_total(counter_total: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_total.launch_counting2(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "d6d0840a27e4782828edde0dcff78f9a"


def test_launch_counting2_smart_shared(counter_smart_shared: Counter, datadir: Path, tmp_path: Path):
    bamfile = datadir / "total.bam"
    countfile = tmp_path / "count.tsv"
    counter_smart_shared.launch_counting2(bamfile, countfile)
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "a09c94e5402521d621a4d94d615cb924"


def test_execute(counter_smart_shared: Counter, tmp_path: Path):
    counter_smart_shared.execute()
    workflow = tmp_path / "workflow.ini"
    assert workflow.exists()
    part1 = tmp_path / "part1/part1.tsv"
    assert part1.exists()
    with part1.open("rb") as out:
        assert md5(out.read()).hexdigest() == "fd9b52548ee66d2c29d7703a20a912bc"
