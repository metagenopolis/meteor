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

"""Test meteor resistance and mapping --kr"""

import json
import random
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from ..counter import Counter
from ..resistance import (
    COVERAGE_BASES,
    SYNONYMOUS_SITES,
    Resistance,
    gene_variants,
    is_synonymous,
)
from ..session import Component


def counts_of(sequence: str, depth: int) -> np.ndarray:
    """Base counts of `depth` reads all carrying the sequence"""
    counts = np.zeros((4, len(sequence)), dtype=np.int64)
    for pos, base in enumerate(sequence):
        counts[COVERAGE_BASES.index(base), pos] = depth
    return counts


def test_synonymous_sites() -> None:
    assert SYNONYMOUS_SITES["TTT"] == pytest.approx(1 / 3)  # TTC
    assert SYNONYMOUS_SITES["GGG"] == pytest.approx(1)  # GGN
    assert SYNONYMOUS_SITES["ATG"] == 0
    assert "TAA" not in SYNONYMOUS_SITES
    assert is_synonymous("AAATTT", 2, "A", "G")  # AAA/AAG: K
    assert not is_synonymous("AAATTT", 3, "T", "C")  # TTT/CTT: F/L
    assert is_synonymous("AAATT", 4, "T", "C") is None  # incomplete codon


def test_gene_variants() -> None:
    sequence = "ATGAAATTTGGGTAA"
    counts = counts_of(sequence, 20)
    counts[:, 0] = [5, 0, 0, 0]  # site below the minimum depth
    counts[:, 5] = [0, 0, 20, 0]  # AAA -> AAG, major synonymous
    counts[:, 6] = [0, 20, 0, 0]  # TTT -> CTT, major non-synonymous
    counts[:, 9] = [5, 0, 15, 0]  # GGG / AGG, minor non-synonymous
    counts[:, 11] = [0, 4, 16, 0]  # GGG / GGC, minor synonymous
    counts[:, 13] = [19, 0, 0, 1]  # 1 read: below the minimum allele count
    result = gene_variants(sequence, counts, min_depth=10, min_alt_count=3, min_frequency=0.1)
    assert result["covered_bp"] == 14
    assert (result["major_snv"], result["major_syn"], result["major_nonsyn"]) == (2, 1, 1)
    assert (result["minor_snv"], result["minor_syn"], result["minor_nonsyn"]) == (2, 1, 1)
    # codons AAA, TTT, GGG (ATG is not covered, TAA is a stop codon)
    assert result["syn_sites"] == pytest.approx(1 / 3 + 1 / 3 + 1)
    assert result["nonsyn_sites"] == pytest.approx(9 - 5 / 3)


def test_gene_variants_out_of_frame() -> None:
    sequence = "ATGAAATTTG"
    counts = counts_of(sequence, 20)
    counts[:, 5] = [0, 0, 20, 0]
    result = gene_variants(sequence, counts, 10, 3, 0.1)
    assert result["major_snv"] == 1
    assert np.isnan(result["major_syn"]) and np.isnan(result["syn_sites"])


# synthetic catalogue: gene 1 = resistance gene of MSP M1, genes 2 and 3 =
# marker genes of M1, gene 4 = resistance gene of no MSP, gene 5 = marker of M2
LENGTH = 300


def _catalogue(
    root: Path, shared: bool = False, duplicate: bool = False, mustard_prediction: bool = True
) -> dict[int, str]:
    rng = random.Random(1)
    genes = {gene: "".join(rng.choice("ACGT") for _ in range(LENGTH)) for gene in range(1, 6)}
    (root / "fasta").mkdir(parents=True)
    (root / "database").mkdir()
    fasta = root / "fasta" / "cat.fasta"
    fasta.write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    pysam.faidx(str(fasta))
    msp_rows = [
        ["M1", 2, "g2", "core"],
        ["M1", 3, "g3", "core"],
        ["M1", 1, "g1", "accessory"],
        ["M2", 5, "g5", "core"],
    ]
    if shared:
        msp_rows.append(["M2", 1, "g1", "accessory"])  # gene 1 belongs to both MSPs
    if duplicate:
        msp_rows.append(msp_rows[1])  # the same marker gene listed twice
    pd.DataFrame(
        msp_rows, columns=["msp_name", "gene_id", "gene_name", "gene_category"]
    ).to_csv(root / "database" / "msp.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"gene_id": [1, 4], "gene_name": ["g1", "g4"], "annotation": ["tet(W)", "blaTEM-1"]}
    ).to_csv(root / "database" / "resfinder.tsv", sep="\t", index=False)
    mustard = pd.DataFrame(
        {"gene_id": [1, 1], "gene_name": ["g1", "g1"], "annotation": ["tetM", "tetM"]}
    )
    if mustard_prediction:
        mustard["Prediction"] = ["Likely_ref", "Likely_ref"]
    mustard.to_csv(root / "database" / "mustard.tsv", sep="\t", index=False)
    # gene 3, marker gene of M1, is gyrA
    pd.DataFrame(
        {"gene_id": [3, 3, 2], "gene_name": ["g3", "g3", "g2"], "annotation": ["K02469", "K00001", "K00002"]}
    ).to_csv(root / "database" / "kegg.tsv", sep="\t", index=False)
    (root / "cat_reference.json").write_text(
        json.dumps(
            {
                "reference_info": {
                    "reference_name": "mock",
                    "database_type": "complete",
                    "reference_date": "2023-03-14",
                },
                "reference_file": {
                    "database_dir": "database",
                    "fasta_dir": "fasta",
                    "fasta_filename": "cat.fasta",
                },
                "annotation": {
                    "msp": {"filename": "msp.tsv"},
                    "resfinder": {"filename": "resfinder.tsv"},
                    "mustard": {"filename": "mustard.tsv"},
                    "kegg": {"filename": "kegg.tsv"},
                },
            }
        )
    )
    return genes


MARKER_GENES = {2, 3, 5}  # core genes of M1 and M2
RESISTANCE_GENES = {1, 3, 4}  # resfinder, mustard and kegg (gyrA) annotations


def _mapped_sample(
    root: Path,
    genes: dict[int, str],
    variants: dict,
    kr: bool = True,
    mapped: tuple = (1, 2, 3, 4),
    split: bool = True,
) -> Path:
    """CRAMs of 100 bp reads every 2 bp on genes 1-4 (depth 50); variants:
    {(gene, pos): (base, every)} gives `base` to one read out of `every`.

    split: as meteor mapping --kr, S1.cram holds the alignments on the marker
    genes and S1_resistance.cram those on the resistance genes (gene 3 in both);
    otherwise S1.cram holds every alignment and the census no resistance file.
    """
    sample = root / "mapping" / "S1"
    sample.mkdir(parents=True)
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": str(g), "LN": LENGTH} for g in genes]}
    )
    files = {"S1": MARKER_GENES, "S1_resistance": RESISTANCE_GENES} if split else {
        "S1": set(genes)
    }
    outputs = {
        name: pysam.AlignmentFile(str(sample / f"{name}.bam"), "wb", header=header)
        for name in files
    }
    n = 0
    for gene in mapped:
        for start in range(0, LENGTH - 100 + 1, 2):
            seq = list(genes[gene][start : start + 100])
            for (vgene, pos), (base, every) in variants.items():
                if vgene == gene and start <= pos < start + 100 and n % every == 0:
                    seq[pos - start] = base
            read = pysam.AlignedSegment(header)
            read.query_name = f"r{n}"
            read.reference_id = gene - 1
            read.reference_start = start
            read.cigarstring = "100M"
            read.query_sequence = "".join(seq)
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.mapping_quality = 1 if gene == 4 else 42
            for name, kept in files.items():
                if gene in kept:
                    outputs[name].write(read)
            n += 1
    for out in outputs.values():
        out.close()
    reference = root / "catalogue" / "fasta" / "cat.fasta"
    for name in files:
        bam = sample / f"{name}.bam"
        cram = sample / f"{name}.cram"
        pysam.view("-C", "-T", str(reference), "-o", str(cram), str(bam), catch_stdout=False)
        pysam.index(str(cram))
        bam.unlink()
    counting = {"resistance_alignments": kr}
    if split and kr:
        counting["resistance_file"] = "S1_resistance.cram"
    (sample / "S1_census_stage_1.json").write_text(
        json.dumps(
            {
                "sample_info": {"sample_name": "S1"},
                "mapping": {"reference_name": "cat"},
                "counting": counting,
            }
        )
    )
    return sample


@pytest.fixture
def resistance(tmp_path: Path) -> tuple[Resistance, dict[int, str]]:
    genes = _catalogue(tmp_path / "catalogue")
    meteor = Component
    meteor.ref_dir = tmp_path / "catalogue"
    meteor.resistance_dir = tmp_path / "resistance"
    meteor.tmp_path = tmp_path / "tmp"
    meteor.threads = 1
    return Resistance(meteor, 10, 3, 0.1, 0.8, 20, 100), genes


def _switch(genes: dict[int, str], gene: int, pos: int) -> str:
    """A base different from the catalogue base"""
    return "A" if genes[gene][pos] != "A" else "C"


def _analysis(root: Path, **catalogue_options) -> tuple[Resistance, dict[int, str]]:
    genes = _catalogue(root, **catalogue_options)
    meteor = Component
    meteor.ref_dir = root
    meteor.resistance_dir = root.parent / "resistance"
    meteor.tmp_path = root.parent / "tmp"
    meteor.threads = 1
    return Resistance(meteor, 10, 3, 0.1, 0.8, 20, 100), genes


def test_gene_variants_tie_favors_the_catalogue_base() -> None:
    sequence = "ATGAAATTTGGGTAA"
    counts = counts_of(sequence, 20)
    counts[:, 3] = [10, 0, 10, 0]  # A / C, an even split: not a major variant
    counts[:, 6] = [0, 15, 0, 5]  # T -> C: a major variant
    result = gene_variants(sequence, counts, 10, 3, 0.1)
    assert result["major_snv"] == 1
    assert result["minor_snv"] == 2  # the even split and the second base of the major site
    assert result["major_nonsyn"] == 1  # AAA -> CAA
    assert result["minor_nonsyn"] == 2  # ATG -> CTG, AAA -> CAA


def test_execute_controls_of_a_multi_msp_gene(tmp_path: Path) -> None:
    """A gene shared by two MSPs sums the controls of both MSPs, without itself"""
    analysis, genes = _analysis(tmp_path / "catalogue", shared=True)
    variants = {(2, 200): (_switch(genes, 2, 200), 1)}  # major variant on a control
    analysis.meteor.mapped_sample_dir = _mapped_sample(
        tmp_path, genes, variants, mapped=(1, 2, 3, 4, 5)
    )
    analysis.execute()
    table = pd.read_table(tmp_path / "resistance" / "S1" / "S1_resistance.tsv")
    arg = table.set_index("gene_id")
    assert arg.loc[1, "msp_name"] == "M1;M2"
    # marker genes of M1 (genes 2 and 3) + marker gene of M2 (gene 5), without gene 1
    assert arg.loc[1, "control_genes"] == 3
    assert arg.loc[1, "control_major_per_kb"] == pytest.approx(
        1000 / arg.loc[1, "control_covered_bp"]
    )
    # gene 3 (gyrA) belongs to M1 only: its controls are the other marker, gene 2
    assert arg.loc[3, "msp_name"] == "M1"
    assert arg.loc[3, "control_genes"] == 1
    # the variant on gene 2, alone on 1/3 of the coverage of the controls of gene 1
    assert arg.loc[3, "control_major_per_kb"] == pytest.approx(3 * arg.loc[1, "control_major_per_kb"])


def test_execute_survives_a_duplicated_msp_gene(tmp_path: Path) -> None:
    """A duplicated (MSP, gene) row in the catalogue must not break the controls"""
    analysis, genes = _analysis(tmp_path / "catalogue", duplicate=True)
    analysis.meteor.mapped_sample_dir = _mapped_sample(tmp_path, genes, {})
    analysis.execute()
    table = pd.read_table(tmp_path / "resistance" / "S1" / "S1_resistance.tsv")
    arg = table.set_index("gene_id")
    # genes 2 and 3, each counted once
    assert arg.loc[1, "control_genes"] == 2
    assert arg.loc[3, "control_genes"] == 1
    assert arg.loc[3, "control_covered_bp"] == arg.loc[1, "control_covered_bp"] - arg.loc[3, "covered_bp"]


def test_execute_without_mustard_prediction_and_zero_read_gene(tmp_path: Path) -> None:
    """A resistance gene without reads is absent from the table but counted in the
    census; a MUSTARD annotation without a Prediction column still runs"""
    analysis, genes = _analysis(tmp_path / "catalogue", mustard_prediction=False)
    analysis.meteor.mapped_sample_dir = _mapped_sample(tmp_path, genes, {}, mapped=(1, 2, 3))
    analysis.execute()
    table = pd.read_table(tmp_path / "resistance" / "S1" / "S1_resistance.tsv")
    arg = table.set_index("gene_id")
    assert table["gene_id"].tolist() == [1, 3]  # gene 4 has no alignment
    assert arg.loc[1, "mustard"] == "tetM" and pd.isna(arg.loc[1, "mustard_prediction"])
    census = json.loads(
        (tmp_path / "resistance" / "S1" / "S1_census_stage_resistance.json").read_text()
    )
    assert census["resistance"]["genes_with_reads"] == 2
    assert census["resistance"]["resistance_genes_total"] == 3


@pytest.mark.parametrize("threads", [1, 2])
def test_execute(resistance, tmp_path: Path, threads: int) -> None:
    analysis, genes = resistance
    analysis.meteor.threads = threads
    variants = {
        (1, 150): (_switch(genes, 1, 150), 1),  # every read: major variant
        (1, 120): (_switch(genes, 1, 120), 4),  # 1 read out of 4: minor variant
        (2, 200): (_switch(genes, 2, 200), 1),  # control: major variant
    }
    analysis.meteor.mapped_sample_dir = _mapped_sample(tmp_path, genes, variants)
    analysis.execute()
    table = pd.read_table(tmp_path / "resistance" / "S1" / "S1_resistance.tsv")
    # genes 1, 3 (gyrA) and 4 have reads, gene 5 (marker of M2) is not a resistance gene
    assert table["gene_id"].tolist() == [1, 3, 4]
    arg = table.set_index("gene_id")
    assert arg.loc[1, "msp_name"] == "M1" and pd.isna(arg.loc[4, "msp_name"])
    assert arg.loc[1, "resfinder"] == "tet(W)" and arg.loc[1, "mustard"] == "tetM"
    assert arg.loc[1, "mustard_prediction"] == "Likely_ref"
    assert arg.loc[4, "resfinder"] == "blaTEM-1" and pd.isna(arg.loc[4, "mustard"])
    assert arg.loc[3, "target"] == "gyrA" and pd.isna(arg.loc[3, "resfinder"])
    assert pd.isna(arg.loc[1, "target"])
    # gene 3 is a marker of M1: its controls are the other marker, gene 2
    assert arg.loc[3, "control_genes"] == 1
    assert arg.loc[3, "control_covered_bp"] == arg.loc[1, "control_covered_bp"] - arg.loc[3, "covered_bp"]
    assert arg.loc[1, "major_snv"] == 1 and arg.loc[1, "minor_snv"] == 1
    assert arg.loc[4, "major_snv"] == 0 and arg.loc[4, "minor_snv"] == 0
    assert arg.loc[1, "major_syn"] + arg.loc[1, "major_nonsyn"] == 1
    # the ends of the genes are covered by fewer than 10 reads
    assert 0.8 < arg.loc[1, "covered_fraction"] < 1
    assert arg.loc[1, "interpretable"]
    assert arg.loc[1, "major_per_kb"] == pytest.approx(1000 / arg.loc[1, "covered_bp"])
    assert arg.loc[1, "multimapped_fraction"] == 0 and arg.loc[4, "multimapped_fraction"] == 1
    # controls: the 2 marker genes of M1, one variant
    assert arg.loc[1, "control_genes"] == 2
    assert arg.loc[1, "control_major_per_kb"] == pytest.approx(
        1000 / arg.loc[1, "control_covered_bp"]
    )
    assert arg.loc[4, "control_genes"] == 0 and pd.isna(arg.loc[4, "control_major_per_kb"])
    census = json.loads(
        (tmp_path / "resistance" / "S1" / "S1_census_stage_resistance.json").read_text()
    )
    assert census["resistance"]["genes_with_reads"] == 3
    assert census["resistance"]["annotations"] == ["resfinder", "mustard", "kegg"]


def test_execute_single_cram_layout(tmp_path: Path) -> None:
    """Resistance genes in <sample>.cram (no resistance file in the census): same
    table as with <sample>_resistance.cram"""
    tables = []
    for split in (True, False):
        root = tmp_path / ("split" if split else "single")
        analysis, genes = _analysis(root / "catalogue")
        variants = {(1, 150): (_switch(genes, 1, 150), 1), (2, 200): (_switch(genes, 2, 200), 1)}
        analysis.meteor.mapped_sample_dir = _mapped_sample(root, genes, variants, split=split)
        analysis.execute()
        tables.append(pd.read_table(root / "resistance" / "S1" / "S1_resistance.tsv"))
    pd.testing.assert_frame_equal(tables[0], tables[1])


def test_streaming_counter_resistance_file(tmp_path: Path) -> None:
    """--kr: each alignment goes to the file of its gene, to both for a gene that
    is a marker gene and a resistance gene; headers list only their genes"""
    from ..counter import StreamingCounter

    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": str(g), "LN": 1000} for g in (1, 2, 3, 4)]}
    )
    streamer = StreamingCounter(
        header, "smart_shared", 0.95,
        tmp_path / "markers.bam", {2, 3},
        tmp_path / "resistance.bam", {1, 3},
    )
    for k, gene in enumerate((1, 2, 3, 4)):
        aln = pysam.AlignedSegment(header)
        aln.query_name = f"r{k}"
        aln.reference_name = str(gene)
        aln.reference_start = 10
        aln.cigarstring = "100M"
        aln.query_sequence = "A" * 100
        aln.set_tag("NM", 0)
        streamer.feed(aln)
    streamer.finish()
    streamer.strain_out.close()
    streamer.resistance_out.close()
    for name, expected in (("markers", ["2", "3"]), ("resistance", ["1", "3"])):
        with pysam.AlignmentFile(str(tmp_path / f"{name}.bam")) as bam:
            assert list(bam.references) == expected
            assert [read.reference_name for read in bam] == expected
    assert streamer.counted_reads == 4


def test_execute_without_kr(resistance, tmp_path: Path) -> None:
    analysis, genes = resistance
    analysis.meteor.mapped_sample_dir = _mapped_sample(tmp_path, genes, {}, kr=False)
    with pytest.raises(SystemExit):
        analysis.execute()


def test_resistance_genes_no_annotation(resistance, tmp_path: Path) -> None:
    analysis, _ = resistance
    ref_json = json.loads((tmp_path / "catalogue" / "cat_reference.json").read_text())
    ref_json["annotation"] = {"msp": ref_json["annotation"]["msp"]}
    with pytest.raises(SystemExit):
        analysis.resistance_genes(ref_json)


# mapping --kr on the test catalogue of test_counter


@pytest.fixture
def counter_kr(tmp_path: Path) -> Counter:
    source = Path(__file__).parent / "test_counter"
    catalogue = tmp_path / "mock"
    shutil.copytree(source / "catalogue" / "mock", catalogue)
    shutil.copytree(source / "part1", tmp_path / "part1")
    meteor = Component
    meteor.tmp_path = tmp_path / "tmp"
    meteor.fastq_dir = tmp_path / "part1"
    meteor.ref_dir = catalogue
    meteor.threads = 1
    meteor.mapping_dir = tmp_path / "mapping"
    meteor.mapping_dir.mkdir()
    return Counter(
        meteor,
        counting_type="smart_shared",
        mapping_type="end-to-end",
        trim=80,
        identity_user=0.95,
        alignment_number=10000,
        core_size=2,
        keep_resistance_alignments=True,
    )


def test_mapping_kr(counter_kr: Counter, monkeypatch) -> None:
    """--kr adds the resistance genes to the filtered alignments of --kf"""
    catalogue = counter_kr.meteor.ref_dir
    ref_json = json.loads((catalogue / "mock_reference.json").read_text())
    markers = counter_kr.get_strain_gene_ids(ref_json)
    # genes with alignments: those of the filtered CRAM of every gene
    sample_dir = counter_kr.meteor.mapping_dir / "part1"
    reference = catalogue / "fasta" / "mock.fasta.gz"
    monkeypatch.setattr(
        Counter, "get_strain_gene_ids", lambda self, ref: set(range(1, 30000))
    )
    counter_kr.keep_resistance_alignments = False
    counter_kr.execute()
    with pysam.AlignmentFile(
        str(sample_dir / "part1.cram"), "rc", reference_filename=str(reference)
    ) as cram:
        aligned = sorted({int(read.reference_name) for read in cram} - markers)
    shutil.rmtree(sample_dir)
    monkeypatch.undo()
    resistance_genes = aligned[:5]
    counter_kr.keep_resistance_alignments = True
    counter_kr.__post_init__()  # new temporary directory
    pd.DataFrame(
        {"gene_id": resistance_genes, "annotation": ["tet(W)"] * len(resistance_genes)}
    ).to_csv(catalogue / "database" / "mock_resfinder.tsv", sep="\t", index=False)
    ref_json["annotation"]["resfinder"] = {"filename": "mock_resfinder.tsv"}
    (catalogue / "mock_reference.json").write_text(json.dumps(ref_json))
    assert counter_kr.keep_filtered_alignments
    counter_kr.execute()
    census = json.loads((sample_dir / "part1_census_stage_1.json").read_text())
    assert census["counting"]["resistance_alignments"] is True
    assert census["counting"]["resistance_file"] == "part1_resistance.cram"
    assert (sample_dir / "part1_resistance.cram.crai").exists()
    found = {}
    for name in ("part1.cram", "part1_resistance.cram"):
        with pysam.AlignmentFile(
            str(sample_dir / name), "rc", reference_filename=str(reference)
        ) as cram:
            found[name] = {int(read.reference_name) for read in cram}
    # part1.cram is the --kf CRAM: marker genes only
    assert found["part1.cram"] and found["part1.cram"] <= markers
    assert found["part1_resistance.cram"] == set(resistance_genes)


def test_mapping_kr_no_annotation(counter_kr: Counter) -> None:
    with pytest.raises(SystemExit):
        counter_kr.execute()
