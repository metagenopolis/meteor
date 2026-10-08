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
import numpy as np
import json

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

def compute_dict_md5(d):
    return md5(json.dumps(d, sort_keys=True).encode("utf-8")).hexdigest()

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
    cram = stage1_dir / "part1_raw.cram"
    assert cram.exists()
    # Compute md5 of alignments in the cram file
    # Cannot compute md5 of the entire file as header may change
    h = md5()
    with AlignmentFile(str(cram.resolve()), "rc") as cramdesc:
        for element in cramdesc:
            h.update(str(element).encode())
    assert h.hexdigest() == "65cdcef09a3f82403765dbbef6e1ae0c"


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
    # Convert pysam alignments to str
    reads.update({k: list(map(str,v)) for k, v in reads.items()})
    assert compute_dict_md5(reads) == "4c46326f906fcebe31bb0197e8882b69"
    assert compute_dict_md5(genes) == "d7e9e2b9b93ad7c449ca90d560de1bd5"


def test_uniq_from_mult(counter_unique: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        reads, genes = counter_unique.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        (unique_reads, genes_mult, unique_on_gene) = counter_unique.uniq_from_mult(
            reads, genes, database
        )
    # Convert pysam alignments to str
    unique_reads.update({k: list(map(str,v)) for k, v in reads.items()})
    assert compute_dict_md5(unique_reads) == "4c46326f906fcebe31bb0197e8882b69"
    assert compute_dict_md5(genes_mult) == "a4161a7477c9610a48669f4834eff776"
    assert compute_dict_md5(unique_on_gene) == "b6d3ba88bdcfcd6dec7cb14d4e01e954"


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
    read_dict.update({k: v.sort() for k, v in read_dict.items()})
    assert compute_dict_md5(read_dict) == "0e6bd9cb3694fb22fa31721513d55fca"
    co = {str(k): v for k, v in co.items()}
    assert compute_dict_md5(co) == "8be9253e6316c77e426ae2ea84940c1c"


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
    # Round as float representation changes between Python <= 3.11 and >= 3.12
    multiple_dict.update({k: round(v,2) for k, v in multiple_dict.items()})
    assert compute_dict_md5(multiple_dict) == "acd524cbfe76351949c071119464a028"


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
        abundance = counter_smart_shared.compute_abs(database, unique_on_gene, multiple_dict)
    # Round as float representation changes between Python <= 3.11 and >= 3.12
    abundance.update({k: round(v,2) for k, v in abundance.items()})
    assert compute_dict_md5(abundance) == "5abede426b36ba12ec85ca24a52a74cc"

def test_compute_abs_total(counter_total: Counter, datadir: Path) -> None:
    cramfile = datadir / "total_raw.cram"
    with AlignmentFile(str(cramfile.resolve()), "rc") as cramdesc:
        _, genes = counter_total.filter_alignments(cramdesc)
        references = map(int, cramdesc.references)
        # get reference length
        lengths = cramdesc.lengths
        database = dict(zip(references, lengths))
        abundance = counter_total.compute_abs_total(database, genes)
    # Round as float representation changes between Python <= 3.11 and >= 3.12
    abundance.update({k: round(v,2) for k, v in abundance.items()})
    assert compute_dict_md5(abundance) == "5089a67644c53346e416ab43ae6d832c"


def test_write_stat(counter_smart_shared: Counter, tmp_path: Path) -> None:
    test = tmp_path / "test.tsv.xz"
    abundance_dict = {11003: 1.5}
    database = {11003: 500}
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
        cram_header = cramdesc.header
    read_list = reads.values()
    merged_list = chain.from_iterable(read_list)
    tmpcramfile = tmp_path / "test"
    counter_unique.save_cram_strain(tmpcramfile, cram_header, merged_list, ref_json)
    assert tmpcramfile.exists()
    # Compute md5 of alignments in the cram file
    # Cannot compute md5 of the entire file as header may change
    h = md5()
    with AlignmentFile(str(tmpcramfile.resolve()), "rc") as cramdesc:
        for element in cramdesc:
            h.update(str(element).encode())
    assert h.hexdigest() == "425664632c57ffe42a93c27f263e0ac5"


def test_launch_counting_unique(counter_unique: Counter, datadir: Path, tmp_path: Path):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = counter_unique.meteor.fastq_dir / "part1_census_stage_1.json"
    census_json = counter_unique.read_json(census_json_file)
    ref_json = counter_unique.read_json(
        counter_unique.meteor.ref_dir / "mock_reference.json"
    )
    counter_unique.launch_counting(
        raw_cramfile, cramfile, countfile, ref_json, census_json, census_json_file
    )
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "f5bc528dcbf594b5089ad7f6228ebab5"
    census_json = counter_unique.read_json(census_json_file)
    assert census_json["counting"]["counted_reads"] == 14386
    assert census_json["counting"]["final_mapping_rate"] == 71.93


def test_launch_counting_total(counter_total: Counter, datadir: Path, tmp_path: Path):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = counter_total.meteor.fastq_dir / "part1_census_stage_1.json"
    census_json = counter_total.read_json(census_json_file)
    ref_json = counter_total.read_json(
        counter_total.meteor.ref_dir / "mock_reference.json"
    )
    counter_total.launch_counting(
        raw_cramfile, cramfile, countfile, ref_json, census_json, census_json_file
    )
    with countfile.open("rb") as out:
        assert md5(out.read()).hexdigest() == "f010e4136323ac408d4c127e243756c2"
    census_json = counter_total.read_json(census_json_file)
    assert census_json["counting"]["counted_reads"] == 14438
    assert census_json["counting"]["final_mapping_rate"] == 72.19


def test_launch_counting_smart_shared(
    counter_smart_shared: Counter, datadir: Path, tmp_path: Path
):
    raw_cramfile = datadir / "total_raw.cram"
    cramfile = datadir / "total.cram"
    countfile = tmp_path / "count.tsv.xz"
    census_json_file = (
        counter_smart_shared.meteor.fastq_dir / "part1_census_stage_1.json"
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
    census_json = counter_smart_shared.read_json(census_json_file)
    assert census_json["counting"]["counted_reads"] == 14438
    assert census_json["counting"]["final_mapping_rate"] == 72.19


def test_execute(counter_smart_shared: Counter, tmp_path: Path):
    counter_smart_shared.execute()
    part1 = tmp_path / "part1/part1.tsv.xz"
    assert part1.exists()
    with part1.open("rb") as out:
        assert md5(out.read()).hexdigest() == "5db950a4404793f73ba034e99cb676fa"


def _toy_alignments(tmp_path: Path) -> Path:
    """Reads of a pair mapped as single-end reads share their name"""
    from pysam import AlignedSegment, AlignmentHeader

    header = AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": str(g), "LN": 1000} for g in (1, 2, 3, 4)]}
    )
    # (name, gene, secondary, mismatches)
    records = [
        ("a", 1, False, 0),  # R1 of pair a: unique on gene 1
        ("a", 1, False, 0),  # R2 of pair a (contiguous): unique on gene 1
        ("b", 2, False, 0),  # R1 of pair b: multi on genes 2 and 3
        ("b", 3, True, 0),
        ("c", 4, False, 0),  # read c: twice on gene 4
        ("c", 4, True, 0),
        ("b", 2, False, 2),  # R2 of pair b (later in the stream): unique on gene 2
        ("b", 3, True, 3),  # lower identity: ignored
    ]
    path = tmp_path / "toy.bam"
    with AlignmentFile(str(path), "wb", header=header) as out:
        for name, gene, secondary, nm in records:
            aln = AlignedSegment(header)
            aln.query_name = name
            aln.reference_name = str(gene)
            aln.reference_start = 10
            aln.flag = 256 if secondary else 0
            aln.cigarstring = "100M"
            aln.query_sequence = "A" * 100
            aln.set_tag("NM", nm)
            out.write(aln)
    return path


def test_mates_sharing_a_name_are_counted_separately(tmp_path: Path) -> None:
    from ..counter import StreamingCounter

    toy = _toy_alignments(tmp_path)
    expected = {
        "unique": ({1: 2, 2: 1, 3: 0, 4: 0}, 3),
        "total": ({1: 2, 2: 2, 3: 1, 4: 2}, 5),
        # b(R1) is split 1/1 on genes 2/3 by unique counts (2: 1, 3: 0);
        # c has no unique evidence: 1/2 per alignment, both on gene 4
        "smart_shared": ({1: 2, 2: 2.0, 3: 0, 4: 1.0}, 5),
    }
    for counting_type, (abundance, counted) in expected.items():
        with AlignmentFile(str(toy)) as cram:
            streamer = StreamingCounter(cram.header, counting_type, 0.95)
            for element in cram:
                streamer.feed(element)
        assert streamer.finish() == abundance, counting_type
        assert streamer.counted_reads == counted, counting_type


def test_filter_alignments_mates_sharing_a_name(
    counter_smart_shared: Counter, tmp_path: Path
) -> None:
    toy = _toy_alignments(tmp_path)
    with AlignmentFile(str(toy)) as cram:
        _, genes = counter_smart_shared.filter_alignments(cram)
    assert sorted(genes.values()) == [[1], [1], [2], [2, 3], [4, 4]]
    database = {1: 1000, 2: 1000, 3: 1000, 4: 1000}
    _, genes_mult, unique_on_gene = counter_smart_shared.uniq_from_mult(
        {k: [] for k in genes}, genes, database
    )
    read_dict, co = counter_smart_shared.compute_co(genes_mult, unique_on_gene)
    multiple = counter_smart_shared.compute_abm(read_dict, co, database)
    assert counter_smart_shared.compute_abs(database, unique_on_gene, multiple) == {
        1: 2, 2: 2.0, 3: 0, 4: 1.0
    }


@pytest.mark.parametrize("counting_type", ["unique", "total", "smart_shared"])
def test_launch_counting_matches_legacy(
    counter_unique: Counter, datadir: Path, tmp_path: Path, counting_type: str
) -> None:
    """Streaming counting == the previous in-memory implementation"""
    import lzma
    counter = counter_unique
    counter.counting_type = counting_type
    ref_json = counter.read_json(counter.meteor.ref_dir / "mock_reference.json")
    census_file = counter.meteor.fastq_dir / "part1_census_stage_1.json"
    results = {}
    for name, launch in (("streaming", counter.launch_counting), ("legacy", counter.launch_counting_legacy)):
        count_file = tmp_path / f"{name}.tsv.xz"
        launch(datadir / "total_raw.cram", tmp_path / f"{name}.cram", count_file, ref_json,
               counter.read_json(census_file), census_file)
        with lzma.open(count_file, "rt") as table:
            results[name] = (pd.read_csv(table, sep="\t"),
                             counter.read_json(census_file)["counting"]["counted_reads"])
        assert (tmp_path / f"{name}.cram").exists()
    streaming, legacy = results["streaming"], results["legacy"]
    assert streaming[1] == legacy[1]
    assert streaming[0][["gene_id", "gene_length"]].equals(legacy[0][["gene_id", "gene_length"]])
    # sums of the same coefficients, possibly in another order (Python < 3.12)
    assert np.allclose(streaming[0]["value"], legacy[0]["value"], rtol=1e-12, atol=0)


def test_launch_counting_without_filtered_alignments(
    counter_smart_shared: Counter, datadir: Path, tmp_path: Path
) -> None:
    counter_smart_shared.keep_filtered_alignments = False
    ref_json = counter_smart_shared.read_json(counter_smart_shared.meteor.ref_dir / "mock_reference.json")
    census_file = counter_smart_shared.meteor.fastq_dir / "part1_census_stage_1.json"
    counter_smart_shared.launch_counting(
        datadir / "total_raw.cram", tmp_path / "strain.cram", tmp_path / "count.tsv.xz",
        ref_json, counter_smart_shared.read_json(census_file), census_file,
    )
    assert (tmp_path / "count.tsv.xz").exists()
    assert not (tmp_path / "strain.cram").exists()


def test_execute_keep_all_then_recount(counter_smart_shared: Counter, tmp_path: Path) -> None:
    """--ka keeps the raw alignments while counting; a second execute counts
    again from them (no mapping) and gives the same table."""
    counter_smart_shared.keep_all_alignments = True
    counter_smart_shared.execute()
    count_file = tmp_path / "part1" / "part1.tsv.xz"
    raw_cram = tmp_path / "part1" / "part1_raw.cram"
    assert raw_cram.exists()
    first = count_file.read_bytes()
    count_file.unlink()
    # new run on the same mapping directory (a new meteor mapping command)
    import dataclasses
    dataclasses.replace(counter_smart_shared, json_data={}).execute()
    assert count_file.read_bytes() == first
    with count_file.open("rb") as out:
        assert md5(out.read()).hexdigest() == "5db950a4404793f73ba034e99cb676fa"


def test_mapping_without_kept_alignments(counter_total: Counter) -> None:
    """Mapping only (no raw alignments kept, no counting)"""
    ref_json = counter_total.read_json(counter_total.meteor.ref_dir / "mock_reference.json")
    census_json_file = counter_total.meteor.fastq_dir / "part1_census_stage_0.json"
    census_json = counter_total.read_json(census_json_file)
    stage1_dir = counter_total.meteor.mapping_dir / census_json["sample_info"]["sample_name"]
    stage1_dir.mkdir(exist_ok=True, parents=True)
    counter_total.json_data[census_json_file] = {
        "census": census_json,
        "directory": stage1_dir,
        "Stage1FileName": stage1_dir / census_json_file.name.replace("stage_0", "stage_1"),
        "reference": ref_json,
    }
    counter_total._build_mapper(keep_raw=False).execute()
    assert counter_total.json_data[census_json_file]["Stage1FileName"].exists()
    assert not (stage1_dir / "part1_raw.cram").exists()
