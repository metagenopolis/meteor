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

"""Single nucleotide variants of the antibiotic resistance genes of a sample,
against the catalogue, with the marker genes of the same species as controls"""

import logging
import lzma
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import get_context
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
from time import perf_counter
from typing import ClassVar

import numpy as np
import pandas as pd
from pysam import AlignmentFile, FastaFile

from meteor.session import Component, Session
from meteor.variantcalling import VariantCalling

# catalogue annotations of acquired antibiotic resistance genes
RESISTANCE_ANNOTATIONS: tuple[str, ...] = ("resfinder", "resfinderfg", "mustard")
# KEGG orthologs of the chromosomal targets of antibiotics, resistant by point
# mutations: quinolones (gyrA, gyrB, parC, parE), rifampicin (rpoB),
# streptomycin (rpsL) and sulfonamides (folP)
TARGET_KOS: dict[str, str] = {
    "K02469": "gyrA",
    "K02470": "gyrB",
    "K02621": "parC",
    "K02622": "parE",
    "K03043": "rpoB",
    "K02950": "rpsL",
    "K00796": "folP",
}
# catalogue annotation of the known resistance positions of the target genes
POSITIONS_ANNOTATION = "resistance_positions"

# standard genetic code, codons in TCAG order
_BASES = "TCAG"
_AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE: dict[str, str] = {
    a + b + c: _AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
# order of the rows of AlignmentFile.count_coverage
COVERAGE_BASES = "ACGT"
_BASE_INDEX = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate(COVERAGE_BASES):
    _BASE_INDEX[ord(_b)] = _i
    _BASE_INDEX[ord(_b.lower())] = _i


def synonymous_sites(codon: str) -> float:
    """Number of synonymous sites of a codon (Nei & Gojobori 1986): for each
    position, the fraction of the 3 other bases giving the same amino acid
    (a change to a stop codon is non-synonymous)"""
    amino_acid = CODON_TABLE[codon]
    sites = 0.0
    for pos in range(3):
        for base in "ACGT":
            if base != codon[pos]:
                mutant = codon[:pos] + base + codon[pos + 1 :]
                sites += (CODON_TABLE[mutant] == amino_acid) / 3.0
    return sites


SYNONYMOUS_SITES: dict[str, float] = {
    codon: synonymous_sites(codon) for codon, aa in CODON_TABLE.items() if aa != "*"
}


def is_synonymous(sequence: str, pos: int, allele_a: str, allele_b: str) -> bool | None:
    """True if the alleles a and b at `pos` code for the same amino acid in the
    codon of the reference sequence (None outside a complete ACGT codon)"""
    start = pos - pos % 3
    codon = sequence[start : start + 3]
    if len(codon) != 3:
        return None
    offset = pos - start
    codon_a = codon[:offset] + allele_a + codon[offset + 1 :]
    codon_b = codon[:offset] + allele_b + codon[offset + 1 :]
    if codon_a not in CODON_TABLE or codon_b not in CODON_TABLE:
        return None
    return CODON_TABLE[codon_a] == CODON_TABLE[codon_b]


def keep_alignment(read) -> bool:
    """Every alignment kept by the mapping step: a multi-mapped read is kept on
    each gene it reaches with its best identity (secondary alignments)"""
    return not (read.is_unmapped or read.is_qcfail)


def gene_variants(
    sequence: str,
    counts: np.ndarray,
    min_depth: int,
    min_alt_count: int,
    min_frequency: float,
    consensus: bool = False,
) -> dict:
    """Variants of one gene from its base counts (4 x length, rows A, C, G, T).

    consensus: also return the consensus sequence of the sample: the major base
    of each covered site, N elsewhere (same coordinates as the catalogue gene).

    A major variant is a covered site whose most frequent base is not the
    catalogue base. A minor variant (polymorphic site) is a covered site whose
    second most frequent base has >= min_alt_count reads and a frequency
    >= min_frequency. Synonymous or not compares the two alleles in the codon
    of the catalogue gene (genes are CDS in frame).
    """
    length = len(sequence)
    ref_index = _BASE_INDEX[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    depth = counts.sum(axis=0)
    covered = (depth >= min_depth) & (ref_index >= 0)
    order = np.argsort(-counts, axis=0, kind="stable")
    positions = np.arange(length)
    top = order[0]
    ref_count = np.where(ref_index >= 0, counts[np.maximum(ref_index, 0), positions], -1)
    # an even split is not a major variant: the catalogue base wins ties
    tie = (top != ref_index) & (counts[top, positions] == ref_count)
    major = np.where(tie, ref_index, top)
    second = np.where(tie, top, order[1])
    major_count = counts[major, positions]
    second_count = counts[second, positions]
    with np.errstate(divide="ignore", invalid="ignore"):
        second_frequency = np.where(depth > 0, second_count / depth, 0.0)
    major_sites = np.flatnonzero(
        covered & (major != ref_index) & (major_count >= min_alt_count)
    )
    minor_sites = np.flatnonzero(
        covered & (second_count >= min_alt_count) & (second_frequency >= min_frequency)
    )
    upper = sequence.upper()
    in_frame = length % 3 == 0
    result: dict[str, float] = {
        "mean_depth": float(depth.mean()) if length else 0.0,
        "covered_bp": int(covered.sum()),
        "major_snv": len(major_sites),
        "minor_snv": len(minor_sites),
        "major_syn": 0,
        "major_nonsyn": 0,
        "minor_syn": 0,
        "minor_nonsyn": 0,
        "syn_sites": np.nan,
        "nonsyn_sites": np.nan,
    }
    if consensus:
        bases = np.frombuffer(COVERAGE_BASES.encode("ascii"), dtype=np.uint8)
        result["consensus"] = (
            np.where(covered, bases[np.maximum(major, 0)], ord("N"))
            .astype(np.uint8)
            .tobytes()
            .decode("ascii")
        )
    if not in_frame:
        for key in ("major_syn", "major_nonsyn", "minor_syn", "minor_nonsyn"):
            result[key] = np.nan
        return result
    for kind, sites, allele_a, allele_b in (
        ("major", major_sites, ref_index, major),
        ("minor", minor_sites, major, second),
    ):
        for pos in sites:
            synonymous = is_synonymous(
                upper, int(pos), COVERAGE_BASES[allele_a[pos]], COVERAGE_BASES[allele_b[pos]]
            )
            if synonymous is not None:
                result[f"{kind}_{'syn' if synonymous else 'nonsyn'}"] += 1
    # sites of the codons fully covered
    codon_covered = covered[: length - length % 3].reshape(-1, 3).all(axis=1)
    syn_sites = 0.0
    codons = 0
    for c in np.flatnonzero(codon_covered):
        codon = upper[3 * c : 3 * c + 3]
        if codon in SYNONYMOUS_SITES:
            syn_sites += SYNONYMOUS_SITES[codon]
            codons += 1
    result["syn_sites"] = syn_sites
    result["nonsyn_sites"] = 3.0 * codons - syn_sites
    return result


def codon_residues(
    bam: AlignmentFile, gene: str, starts: list[int], base_quality: int
) -> dict[int, dict[str, int]]:
    """Amino acids carried by the reads at the codons starting at `starts`
    (0-based) of a gene, read by read: a codon counts for a read when its 3
    bases are aligned without indel and have a quality >= base_quality.
    Counting each base separately would mix the bases of different reads into
    codons that no read carries.

    :return: {codon start: {amino acid: reads}}
    """
    residues: dict[int, dict[str, int]] = {start: {} for start in starts}
    if not starts:
        return residues
    low, high = min(starts), max(starts) + 3
    for read in bam.fetch(gene, low, high):
        sequence = read.query_sequence
        qualities = read.query_qualities
        if not keep_alignment(read) or sequence is None or qualities is None:
            continue
        query_of = {
            ref: query
            for query, ref in read.get_aligned_pairs(matches_only=True)
            if low <= ref < high
        }
        for start in starts:
            query = query_of.get(start)
            if (
                query is None
                or query_of.get(start + 1) != query + 1
                or query_of.get(start + 2) != query + 2
                or min(qualities[query : query + 3]) < base_quality
            ):
                continue
            amino_acid = CODON_TABLE.get(sequence[query : query + 3].upper())
            if amino_acid is not None:
                counts = residues[start]
                counts[amino_acid] = counts.get(amino_acid, 0) + 1
    return residues


def resistance_worker(
    bam_path: str,
    ref_path: str,
    genes: list[tuple[str, int]],
    min_depth: int,
    min_alt_count: int,
    min_frequency: float,
    base_quality: int,
    consensus: bool = False,
    positions: dict[str, list[int]] | None = None,
) -> list[dict]:
    """Variants of a batch of genes (module level: run in spawn workers);
    positions: {gene: codon starts} whose residues are counted (`residues`)"""
    rows = []
    positions = positions or {}
    with AlignmentFile(bam_path, "rb") as bam, FastaFile(ref_path) as fasta:
        for gene, length in genes:
            alignments = multimapped = 0
            for read in bam.fetch(gene):
                if keep_alignment(read):
                    alignments += 1
                    # bowtie2 gives MAPQ 0 or 1 to reads with equally good hits
                    multimapped += read.mapping_quality <= 1
            counts = np.array(
                bam.count_coverage(
                    gene,
                    0,
                    length,
                    quality_threshold=base_quality,
                    read_callback=keep_alignment,
                ),
                dtype=np.int64,
            )
            row: dict = {"gene_id": int(gene), "gene_length": length, "alignments": alignments}
            row["multimapped_fraction"] = multimapped / alignments if alignments else np.nan
            row.update(
                gene_variants(
                    fasta.fetch(gene), counts, min_depth, min_alt_count, min_frequency, consensus
                )
            )
            if gene in positions:
                row["residues"] = codon_residues(bam, gene, positions[gene], base_quality)
            rows.append(row)
    return rows


def load_kegg_targets(path: Path) -> pd.DataFrame:
    """The genes of the target KEGG orthologs, without loading the whole table:
    the KEGG annotation of a complete catalogue holds millions of rows"""
    kos = set(TARGET_KOS)
    if path.suffix == ".feather":
        import pyarrow.feather as feather

        table = feather.read_feather(path, columns=["gene_id", "annotation"])
        return table[table["annotation"].isin(kos)]
    chunks = pd.read_csv(path, sep="\t", usecols=["gene_id", "annotation"], chunksize=2_000_000)
    return pd.concat((chunk[chunk["annotation"].isin(kos)] for chunk in chunks), ignore_index=True)


def resistance_gene_table(session: Session, ref_dir: Path, ref_json: dict) -> pd.DataFrame:
    """One row per resistance gene of the catalogue: its ResFinder, ResFinderFG
    and MUSTARD annotations, and `target` for the genes of the KEGG orthologs of
    TARGET_KOS (empty table without any of these annotations)"""
    database_dir = ref_dir / ref_json["reference_file"]["database_dir"]
    annotations = ref_json.get("annotation", {})
    columns = []
    for name in RESISTANCE_ANNOTATIONS:
        if name not in annotations:
            continue
        table = session.load_data(database_dir / annotations[name]["filename"])
        table = table.assign(gene_id=table["gene_id"].astype(int))
        columns.append(
            table.groupby("gene_id")["annotation"]
            .agg(lambda values: ",".join(sorted(set(map(str, values)))))
            .rename(name)
        )
        if name == "mustard" and "Prediction" in table:
            columns.append(
                table.groupby("gene_id")["Prediction"]
                .agg(lambda values: ",".join(sorted(set(map(str, values)))))
                .rename("mustard_prediction")
            )
    if "kegg" in annotations:
        kegg = load_kegg_targets(database_dir / annotations["kegg"]["filename"])
        columns.append(
            kegg.assign(
                gene_id=kegg["gene_id"].astype(int),
                target=kegg["annotation"].map(TARGET_KOS),
            )
            .groupby("gene_id")["target"]
            .agg(lambda values: ",".join(sorted(set(values))))
        )
    names = [*RESISTANCE_ANNOTATIONS, "mustard_prediction", "target"]
    if not columns:
        return pd.DataFrame(columns=["gene_id", *names])
    genes = pd.concat(columns, axis=1)
    for name in names:
        if name not in genes:
            genes[name] = pd.NA
    return genes[names].rename_axis("gene_id").reset_index()


def _ratio(numerator, denominator):
    """numerator / denominator, NaN when the denominator is 0"""
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denominator > 0, numerator / denominator, np.nan)


# columns of the catalogue table of the known resistance positions
POSITION_COLUMNS = (
    "gene_id",
    "family",
    "codon_pos",
    "ref_residue",
    "resistance_residues",
    "gene_aa_pos",
    "gene_nt_pos",
    "gene_residue",
    "window_identity",
)


def position_table(
    positions: pd.DataFrame,
    residues: dict[int, dict[int, dict[str, int]]],
    min_depth: int,
    min_alt_count: int,
    min_frequency: float,
) -> pd.DataFrame:
    """Residues of the sample at the known resistance positions of its genes.

    positions: catalogue positions (gene_nt_pos: first base of the codon,
    1-based); residues: {gene_id: {codon start (0-based): {amino acid: reads}}}.
    `known_resistance`: >= min_alt_count reads and >= min_frequency of them
    carry a resistance residue (a minor resistant population is kept), at a
    position read by >= min_depth reads.
    """
    rows = []
    for row in positions.itertuples(index=False):
        counts = residues.get(int(row.gene_id), {}).get(int(row.gene_nt_pos) - 1, {})
        depth = sum(counts.values())
        if depth == 0:
            continue
        # most reads first, then alphabetical: a stable major residue
        ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        resistant = set(str(row.resistance_residues).split(","))
        resistance_reads = sum(n for residue, n in counts.items() if residue in resistant)
        rows.append(
            {
                **{column: getattr(row, column) for column in POSITION_COLUMNS},
                "depth": depth,
                "readable": depth >= min_depth,
                "major_residue": ordered[0][0],
                "major_frequency": ordered[0][1] / depth,
                "major_is_resistance": ordered[0][0] in resistant,
                "resistance_reads": resistance_reads,
                "resistance_frequency": resistance_reads / depth,
                "known_resistance": depth >= min_depth
                and resistance_reads >= min_alt_count
                and resistance_reads / depth >= min_frequency,
                "residue_counts": ",".join(f"{residue}:{n}" for residue, n in ordered),
            }
        )
    columns = [
        *POSITION_COLUMNS, "depth", "readable", "major_residue", "major_frequency",
        "major_is_resistance", "resistance_reads", "resistance_frequency",
        "known_resistance", "residue_counts",
    ]
    return pd.DataFrame(rows, columns=columns)


@dataclass
class Resistance(Session):
    """Variants of the antibiotic resistance genes of a sample (meteor resistance)"""

    DEFAULT_MIN_DEPTH: ClassVar[int] = 10
    DEFAULT_MIN_ALT_COUNT: ClassVar[int] = 3
    DEFAULT_MIN_FREQUENCY: ClassVar[float] = 0.1
    DEFAULT_MIN_GENE_COVERAGE: ClassVar[float] = 0.8
    DEFAULT_BASE_QUALITY: ClassVar[int] = 20
    DEFAULT_NUM_THREADS: ClassVar[int] = 1

    meteor: type[Component]
    min_depth: int
    min_alt_count: int
    min_frequency: float
    min_gene_coverage: float
    base_quality: int
    core_size: int

    def __post_init__(self) -> None:
        if self.meteor.tmp_path:
            self.meteor.tmp_path.mkdir(exist_ok=True, parents=True)
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.resistance_dir.mkdir(exist_ok=True, parents=True)

    def resistance_genes(self, ref_json: dict) -> pd.DataFrame:
        """One row per resistance gene: its annotation in each database"""
        genes = resistance_gene_table(self, self.meteor.ref_dir, ref_json)
        if genes.empty:
            logging.error(
                "The catalogue has no resistance annotation (%s or kegg)",
                ", ".join(RESISTANCE_ANNOTATIONS),
            )
            sys.exit(1)
        return genes

    def gene_stats(
        self,
        bam: Path,
        reference: str,
        genes: list[tuple[str, int]],
        consensus: bool = False,
        positions: dict[str, list[int]] | None = None,
    ) -> tuple[pd.DataFrame, dict[int, dict[int, dict[str, int]]]]:
        """Variants of every gene, in worker processes when threads > 1, and
        the residues of the reads at the codons of `positions` ({gene: starts})

        :return: (one row per gene, {gene_id: {codon start: {amino acid: reads}}})
        """
        positions = positions or {}
        args = (
            self.min_depth, self.min_alt_count, self.min_frequency, self.base_quality, consensus
        )
        workers = max(1, self.meteor.threads)
        if workers == 1 or len(genes) < 2 * workers:
            rows = resistance_worker(str(bam), reference, genes, *args, positions)
        else:
            batch = max(1, min(2000, len(genes) // (workers * 8)))
            batches = [genes[i : i + batch] for i in range(0, len(genes), batch)]
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) as pool:
                results = pool.map(
                    resistance_worker,
                    *zip(
                        *[
                            (
                                str(bam), reference, b, *args,
                                {g: positions[g] for g, _ in b if g in positions},
                            )
                            for b in batches
                        ]
                    ),
                )
                rows = [row for result in results for row in result]
        residues = {row["gene_id"]: row.pop("residues") for row in rows if "residues" in row}
        return pd.DataFrame(rows), residues

    def add_rates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Coverage, rates per kb covered, dN/dS against the catalogue (major
        variants) and pN/pS within the sample (minor variants)"""
        df["covered_fraction"] = _ratio(df["covered_bp"], df["gene_length"])
        df["interpretable"] = df["covered_fraction"] >= self.min_gene_coverage
        df["major_per_kb"] = 1000 * _ratio(df["major_snv"], df["covered_bp"])
        df["minor_per_kb"] = 1000 * _ratio(df["minor_snv"], df["covered_bp"])
        for kind, name in (("major", "major_dnds"), ("minor", "minor_pnps")):
            df[name] = _ratio(
                _ratio(df[f"{kind}_nonsyn"], df["nonsyn_sites"]),
                _ratio(df[f"{kind}_syn"], df["syn_sites"]),
            )
        return df

    def load_positions(self, ref_json: dict, gene_ids: set[int]) -> dict:
        """Known resistance positions of the resistance genes, when the catalogue
        has them: {"table": catalogue rows, "starts": {gene_id: codon starts (0-based)}}"""
        annotation = ref_json.get("annotation", {}).get(POSITIONS_ANNOTATION)
        if annotation is None:
            return {"table": None, "starts": {}}
        database_dir = self.meteor.ref_dir / ref_json["reference_file"]["database_dir"]
        table = self.load_data(database_dir / annotation["filename"])
        table = table.assign(gene_id=table["gene_id"].astype(int))
        table = table[table["gene_id"].isin(gene_ids)][list(POSITION_COLUMNS)]
        starts: dict[int, set[int]] = {}
        for gene, nt_pos in zip(table["gene_id"], table["gene_nt_pos"]):
            starts.setdefault(int(gene), set()).add(int(nt_pos) - 1)
        logging.info("%d known resistance positions on %d genes", len(table), len(starts))
        return {"table": table, "starts": {gene: sorted(s) for gene, s in starts.items()}}

    def execute(self) -> None:
        """Variants of the resistance genes of the mapped sample"""
        try:
            self._execute()
        finally:
            rmtree(self.meteor.tmp_dir, ignore_errors=True)

    def _execute(self) -> None:
        start = perf_counter()
        ref_json = self.get_reference_info(self.meteor.ref_dir)
        Component.check_catalogue(ref_json)
        census = self.get_census_stage(self.meteor.mapped_sample_dir, 1)
        sample = census["sample_info"]["sample_name"]
        if not census.get("counting", {}).get("resistance_alignments", False):
            logging.error(
                "%s was not mapped with --kr: its filtered alignments do not hold the "
                "resistance genes. Map the sample again with meteor mapping --kr.",
                sample,
            )
            sys.exit(1)
        # alignments on the marker genes (controls) and on the resistance genes
        cram = self.meteor.mapped_sample_dir / f"{sample}.cram"
        resistance_cram = self.meteor.mapped_sample_dir / census["counting"].get(
            "resistance_file", cram.name
        )
        reference = (
            self.meteor.ref_dir
            / ref_json["reference_file"]["fasta_dir"]
            / ref_json["reference_file"]["fasta_filename"]
        )
        database_dir = self.meteor.ref_dir / ref_json["reference_file"]["database_dir"]
        resistance = self.resistance_genes(ref_json)
        msp = self.load_data(database_dir / ref_json["annotation"]["msp"]["filename"])
        msp = msp.assign(gene_id=msp["gene_id"].astype(int))
        markers = msp[msp["gene_category"] == "core"].groupby("msp_name").head(self.core_size)
        msp_of_gene = (
            msp[msp["gene_id"].isin(set(resistance["gene_id"]))]
            .groupby("gene_id")["msp_name"]
            .agg(lambda names: ";".join(sorted(set(names))))
        )
        carrying = {name for names in msp_of_gene for name in names.split(";")}
        controls = markers[markers["msp_name"].isin(carrying)]
        resistance_ids = set(resistance["gene_id"])
        positions = self.load_positions(ref_json, resistance_ids)
        helper = VariantCalling(self.meteor, {}, 0, 0, 0, 0.0, 1, self.core_size)
        # each resistance gene is read from the resistance alignments, or from
        # the marker alignments when the resistance CRAM does not list it (a
        # target added to --kr after the mapping: a marker gene is still there);
        # a resistance gene that is also a marker gene is read once
        if resistance_cram == cram:
            # single-file layout: every catalogue gene is in the header
            in_resistance_cram = resistance_ids
        else:
            with AlignmentFile(
                str(resistance_cram.resolve()), "rc", reference_filename=str(reference.resolve())
            ) as alignments:
                in_resistance_cram = resistance_ids & {int(name) for name in alignments.references}
        from_markers = resistance_ids - in_resistance_cram
        parts = []
        residues: dict[int, dict[int, dict[str, int]]] = {}
        for source, targets, with_consensus in (
            (resistance_cram, in_resistance_cram, True),
            (cram, (set(controls["gene_id"]) | from_markers) - in_resistance_cram, bool(from_markers)),
        ):
            if not targets:
                continue
            bam, header_genes = helper.write_marker_alignments(
                source, reference, targets, self.meteor.tmp_dir
            )
            fasta = helper.write_marker_reference(reference, header_genes, self.meteor.tmp_dir)
            with AlignmentFile(str(bam), "rb") as alignments:
                genes = [
                    (name, length)
                    for name, length in zip(alignments.references, alignments.lengths)
                    if int(name) in targets
                ]
            logging.info("Variants of %d genes of %s", len(genes), source.name)
            if genes:
                stats_part, residues_part = self.gene_stats(
                    bam,
                    fasta,
                    genes,
                    consensus=with_consensus,
                    positions={
                        name: starts
                        for name, _ in genes
                        if (starts := positions["starts"].get(int(name)))
                    },
                )
                parts.append(stats_part)
                residues.update(residues_part)
        if not parts:
            logging.error("%s: no alignment on the resistance genes or their controls", sample)
            sys.exit(1)
        stats = self.add_rates(pd.concat(parts, ignore_index=True))
        if "consensus" not in stats:
            stats["consensus"] = None
        consensus = stats.set_index("gene_id")["consensus"]
        stats = stats.drop(columns="consensus")
        # controls: interpretable marker genes of the MSPs carrying the gene
        # (without the gene itself: gyrA is often a marker gene)
        sums = ["covered_bp", "major_snv", "minor_snv"]
        marker_stats = controls.merge(stats, on="gene_id")
        # a duplicated (MSP, gene) row in the catalogue would be counted twice
        # and break the Series lookups below
        marker_stats = marker_stats.drop_duplicates(["msp_name", "gene_id"])
        marker_stats = marker_stats[marker_stats["interpretable"]].assign(control_genes=1)
        per_msp = marker_stats.groupby("msp_name")[[*sums, "control_genes"]].sum()
        own = marker_stats.set_index(["msp_name", "gene_id"])[[*sums, "control_genes"]]
        table = resistance.merge(stats, on="gene_id")
        table = table[table["alignments"] > 0].copy()
        table.insert(1, "msp_name", table["gene_id"].map(msp_of_gene))
        control_rows = []
        for gene, names in zip(table["gene_id"], table["msp_name"]):
            rows = pd.Series(0.0, index=per_msp.columns)
            if isinstance(names, str):
                for name in names.split(";"):
                    if name in per_msp.index:
                        rows = rows + per_msp.loc[name]
                    if (name, gene) in own.index:
                        rows = rows - own.loc[(name, gene)]
            control_rows.append(rows)
        control = pd.DataFrame(
            control_rows,
            index=table.index,
            columns=["covered_bp", "major_snv", "minor_snv", "control_genes"],
        )
        table["control_genes"] = control["control_genes"].astype(int)
        table["control_covered_bp"] = control["covered_bp"].astype(int)
        table["control_major_per_kb"] = 1000 * _ratio(control["major_snv"], control["covered_bp"])
        table["control_minor_per_kb"] = 1000 * _ratio(control["minor_snv"], control["covered_bp"])
        sample_dir = self.meteor.resistance_dir / sample
        sample_dir.mkdir(exist_ok=True, parents=True)
        output = sample_dir / f"{sample}_resistance.tsv"
        table = table.sort_values("gene_id")
        table.to_csv(output, sep="\t", index=False, na_rep="NA")
        # consensus sequences of the interpretable genes (meteor resistance_summary --trees)
        consensus_file = sample_dir / f"{sample}_resistance_consensus.fasta.xz"
        with lzma.open(consensus_file, "wt", encoding="UTF-8") as fasta_out:
            for gene in table.loc[table["interpretable"], "gene_id"]:
                sequence = consensus.get(gene)
                if isinstance(sequence, str):
                    fasta_out.write(f">{gene}\n{sequence}\n")
        # residues at the known resistance positions of the target genes
        positions_file = None
        if positions["table"] is not None:
            position_rows = position_table(
                positions["table"],
                residues,
                self.min_depth,
                self.min_alt_count,
                self.min_frequency,
            )
            position_rows.insert(1, "msp_name", position_rows["gene_id"].map(msp_of_gene))
            positions_file = sample_dir / f"{sample}_resistance_positions.tsv"
            position_rows.sort_values(["gene_id", "codon_pos"]).to_csv(
                positions_file, sep="\t", index=False, na_rep="NA"
            )
            logging.info(
                "%d known resistance positions read, %d with a resistance residue: %s",
                int(position_rows["readable"].sum()),
                int(position_rows["known_resistance"].sum()),
                positions_file,
            )
        annotations = ref_json.get("annotation", {})
        self.save_config(
            {
                "meteor_version": self.meteor.version,
                "sample_info": census["sample_info"],
                "mapping": {"reference_name": ref_json["reference_info"]["reference_name"]},
                "resistance": {
                    "resistance_date": datetime.now().strftime("%Y-%m-%d"),
                    "annotations": [
                        n
                        for n in (*RESISTANCE_ANNOTATIONS, "kegg", POSITIONS_ANNOTATION)
                        if n in annotations
                    ],
                    "target_kos": sorted(TARGET_KOS),
                    "min_depth": self.min_depth,
                    "min_alt_count": self.min_alt_count,
                    "min_frequency": self.min_frequency,
                    "min_gene_coverage": self.min_gene_coverage,
                    "base_quality": self.base_quality,
                    "core_size": self.core_size,
                    "resistance_genes_total": len(resistance),
                    "genes_with_reads": len(table),
                    "interpretable_genes": int(table["interpretable"].sum()),
                    # resistance genes missing from the resistance CRAM, read from
                    # the marker alignments (mapped before they were added to --kr)
                    "genes_read_from_markers": int(table["gene_id"].isin(from_markers).sum()),
                    "resistance_file": output.name,
                    "consensus_file": consensus_file.name,
                    **({"positions_file": positions_file.name} if positions_file else {}),
                },
            },
            sample_dir / f"{sample}_census_stage_resistance.json",
        )
        logging.info(
            "%d resistance genes with reads, %d interpretable: %s",
            len(table),
            int(table["interpretable"].sum()),
            output,
        )
        logging.info("Completed resistance analysis in %f seconds", perf_counter() - start)
