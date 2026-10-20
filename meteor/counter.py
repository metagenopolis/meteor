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

"""Run mapping and performs counting"""

import logging
import sys
import pysam
import lzma
from dataclasses import dataclass, field
from tempfile import mkdtemp, mkstemp
from pathlib import Path
from meteor.mapper import Mapper
from meteor.session import Session, Component
from typing import ClassVar
from collections import defaultdict, Counter as TupleCounter
from itertools import chain, repeat
from pysam import index, sort, AlignmentFile, AlignmentHeader, AlignedSegment  # type: ignore[attr-defined]
from time import perf_counter
from shutil import rmtree


class StreamingCounter:
    """Identity filtering + counting fed one alignment at a time.

    Only the alignments of a read with the best identity (above the identity
    threshold) are kept. "unique" counts the reads reaching one gene, "total"
    counts every kept alignment, "smart_shared" adds to the unique counts each
    multi-mapped read split between its genes in proportion to their unique
    counts (1 / number of alignments when none of them has unique reads).
    Only the AlignedSegment of the current read is kept in memory.

    Every read is counted on its own (single-end counting). bowtie2 writes the
    alignments of a read contiguously and the first one is primary (the -k
    extra ones are flagged secondary), so a read starts at each primary
    alignment or at each change of name. The read name is never used as a key:
    R1 and R2 of a pair mapped as single-end reads share their name (bowtie2
    drops the Illumina comment "1:N:0:..." / "2:N:0:...") but are two reads.
    A read is counted as soon as its last alignment is seen; multi-mapped reads
    are aggregated by gene list, so smart_shared coefficients are computed once
    per distinct list.

    A coordinate-sorted input (raw CRAM sorted by an older pipeline) no longer
    keeps the alignments of a read together: alignments are then grouped by
    name, as in older meteor versions, and reads sharing a name are merged.
    """

    def __init__(
        self,
        header: AlignmentHeader,
        counting_type: str,
        identity_threshold: float,
        strain_path: Path | None = None,
        strain_gene_ids: set[int] | None = None,
    ) -> None:
        self.counting_type = counting_type
        self.identity_threshold = identity_threshold
        references = header.references
        lengths = header.lengths
        self.gene_of_tid = [int(ref) for ref in references]
        self.database: dict[int, int] = dict(zip(self.gene_of_tid, lengths))
        # bytearray indexed by reference id: 1 when the gene is a marker gene
        self.strain_tid = bytearray(len(references))
        self.strain_out: AlignmentFile | None = None
        self.strain_new_tid: list[int] = []
        if strain_path is not None and strain_gene_ids is not None:
            # The filtered alignments get a header restricted to marker genes:
            # with the bowtie2 header (one @SQ per catalogue gene) every later
            # reader/writer handles millions of sequences (and writing CRAM
            # computes the MD5 of every catalogue sequence).
            marker_tids = [
                tid
                for tid, gene in enumerate(self.gene_of_tid)
                if gene in strain_gene_ids
            ]
            self.strain_new_tid = [-1] * len(references)
            for new, tid in enumerate(marker_tids):
                self.strain_tid[tid] = 1
                self.strain_new_tid[tid] = new
            strain_header = AlignmentHeader.from_dict(
                {
                    "HD": {"VN": "1.6", "SO": "unsorted"},
                    "SQ": [
                        {"SN": references[tid], "LN": lengths[tid]}
                        for tid in marker_tids
                    ],
                }
            )
            self.strain_out = AlignmentFile(
                str(strain_path.resolve()), "wb0", header=strain_header
            )
        self.unique_on_gene: dict[int, int] = dict.fromkeys(self.database, 0)
        self.total_on_gene: dict[int, int] = (
            dict.fromkeys(self.database, 0) if counting_type == "total" else {}
        )
        self.multi: TupleCounter = TupleCounter()
        self.counted_reads = 0
        # current read: name, best identity, reference ids reaching it and
        # alignments on marker genes reaching it
        self._read: str | None = None
        self._best = -1.0
        self._tids: list[int] = []
        self._alns: list[AlignedSegment] = []
        # coordinate-sorted input: read name -> [best identity, reference ids,
        # marker alignments]
        self.by_name = header.to_dict().get("HD", {}).get("SO") == "coordinate"
        self.states: dict[str, list] = {}
        if self.by_name:
            logging.warning(
                "Coordinate-sorted alignments: alignments sharing a read name "
                "are merged (R1 and R2 with the same name count as one read). "
                "Map the sample again to count every read."
            )
            self.feed = self._feed_by_name  # type: ignore[method-assign]

    def feed(self, element: AlignedSegment) -> None:
        read_id = element.query_name
        if not element.is_secondary or read_id != self._read:
            # a new read starts (even when its name was already seen)
            if self._tids:
                self._count()
            self._read = read_id
            self._best = -1.0
            self._tids = []
            self._alns = []
        stats = element.get_cigar_stats()[0]
        ali = stats[0] + stats[1] + stats[2]
        identity = (ali - element.get_tag("NM")) / ali
        if identity < self.identity_threshold:
            return
        tid = element.reference_id
        if identity == self._best:
            self._tids.append(tid)
            if self.strain_tid[tid]:
                self._alns.append(element)
        elif identity > self._best:
            self._best = identity
            self._tids = [tid]
            self._alns = [element] if self.strain_tid[tid] else []

    def _feed_by_name(self, element: AlignedSegment) -> None:
        stats = element.get_cigar_stats()[0]
        ali = stats[0] + stats[1] + stats[2]
        identity = (ali - element.get_tag("NM")) / ali
        if identity < self.identity_threshold:
            return
        tid = element.reference_id
        marker = [element] if self.strain_tid[tid] else []
        state = self.states.get(element.query_name)
        if state is None:
            self.states[element.query_name] = [identity, [tid], marker]
        elif identity == state[0]:
            state[1].append(tid)
            state[2].extend(marker)
        elif identity > state[0]:
            state[0] = identity
            state[1] = [tid]
            state[2] = marker

    def _count(self) -> None:
        """Count the current read and write its marker alignments"""
        self._count_read(self._tids, self._alns)

    def _count_read(self, tids: list[int], alns: list[AlignedSegment]) -> None:
        """Count one read and write its marker alignments"""
        gene_of_tid = self.gene_of_tid
        if len(tids) == 1:
            self.unique_on_gene[gene_of_tid[tids[0]]] += 1
        else:
            if self.counting_type == "unique":
                return
            if self.counting_type == "smart_shared":
                self.multi[tuple(gene_of_tid[t] for t in tids)] += 1
        if self.counting_type == "total":
            total_on_gene = self.total_on_gene
            for t in tids:
                total_on_gene[gene_of_tid[t]] += 1
        self.counted_reads += 1
        if self.strain_out is not None and alns:
            new_tid = self.strain_new_tid
            write = self.strain_out.write
            for aln in alns:
                aln.reference_id = new_tid[aln.reference_id]
                if aln.next_reference_id >= 0:
                    aln.next_reference_id = new_tid[aln.next_reference_id]
                write(aln)

    def finish(self) -> dict[int, int | float]:
        """Count the last read and return the abundance of every gene."""
        if self._tids:
            self._count()
        states = self.states
        self.states = {}
        for _, tids, alns in states.values():
            self._count_read(tids, alns)
        del states
        self._read = None
        self._tids = []
        self._alns = []
        if self.counting_type == "total":
            return self.total_on_gene
        if self.counting_type == "unique":
            return dict(self.unique_on_gene)
        return self._smart_shared()

    def _smart_shared(self) -> dict[int, int | float]:
        """Abundance of every gene in smart_shared mode.

        Each read contributes its coefficient once per gene and the contributions
        are added with the built-in sum() (Python >= 3.12 sum() is compensated,
        so the result does not depend on read order).
        """
        unique_on_gene = self.unique_on_gene
        # gene -> {coefficient: number of reads}
        contrib: dict[int, dict[float, int]] = defaultdict(lambda: defaultdict(int))
        for genes, n_reads in self.multi.items():
            som = sum(unique_on_gene[gene] for gene in genes)
            per_gene: dict[int, float] = {}
            if som == 0:
                # 1 / nb of alignments for each alignment (a gene reached
                # twice by the read gets 2 / nb)
                coef = 1.0 / len(genes)
                for gene in genes:
                    per_gene[gene] = per_gene.get(gene, 0.0) + coef
            else:
                seen: set[int] = set()
                duplicated = {gene for gene in genes if gene in seen or seen.add(gene)}
                for gene in genes:
                    nb_unique = unique_on_gene[gene]
                    if nb_unique == 0:
                        continue
                    if gene in duplicated:
                        per_gene[gene] = per_gene.get(gene, 0.0) + nb_unique / float(som)
                    else:
                        per_gene[gene] = nb_unique / float(som)
            for gene, coef in per_gene.items():
                contrib[gene][coef] += n_reads
        abundance: dict[int, int | float] = dict(unique_on_gene)
        for gene, coefs in contrib.items():
            abundance[gene] = unique_on_gene[gene] + sum(
                chain.from_iterable(repeat(c, n) for c, n in coefs.items())
            )
        return abundance


@dataclass
class Counter(Session):
    """Counter session map and count"""

    COUNTING_TYPES: ClassVar[list[str]] = ["total", "smart_shared", "unique"]
    DEFAULT_COUNTING_TYPE: ClassVar[str] = "smart_shared"
    NO_IDENTITY_THRESHOLD: ClassVar[float] = 0.0
    DEFAULT_IDENTITY_THRESHOLD_COMPLETE: ClassVar[float] = 0.95
    DEFAULT_IDENTITY_THRESHOLD_TAXO: ClassVar[float] = 0.98

    meteor: type[Component]
    counting_type: str
    mapping_type: str
    trim: int
    identity_user: float | None
    alignment_number: int
    core_size: int
    keep_all_alignments: bool = False
    keep_filtered_alignments: bool = False
    json_data: dict = field(default_factory=dict)
    identity_threshold: float = field(default_factory=float)

    def __post_init__(self) -> None:
        if self.counting_type not in Counter.COUNTING_TYPES:
            raise ValueError(f"{self.counting_type} is not a valid counting type")

        if self.meteor.tmp_path:
            self.meteor.tmp_path.mkdir(exist_ok=True)
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.mapping_dir.mkdir(exist_ok=True)

    def launch_mapping(self) -> None:
        """Create temporary indexed files and map against"""
        self._build_mapper().execute()

    def _build_mapper(self, keep_raw: bool = True, consumer_factory=None) -> Mapper:
        """Mapper for all the libraries (fastq files) of the sample"""
        fastq_paths = []
        # loop on each library
        for dict_data in self.json_data.values():
            census = dict_data["census"]
            sample_file = census["sample_file"]
            # reindexing this library reads and fill FLibraryIndexerReport
            fastq_paths.append(str(self.meteor.fastq_dir / sample_file["fastq_file"]))

        # mapping this library on the reference
        return Mapper(
            self.meteor,
            dict_data,
            fastq_paths,
            self.mapping_type,
            self.trim,
            self.alignment_number,
            keep_raw,
            consumer_factory,
        )

    def get_strain_gene_ids(self, ref_json: dict) -> set[int]:
        """Core genes (core_size per MSP) kept in the filtered CRAM"""
        msp_file = (
            self.meteor.ref_dir
            / ref_json["reference_file"]["database_dir"]
            / ref_json["annotation"]["msp"]["filename"]
        )
        msp_content = (
            self.load_data(msp_file)
            .query("gene_category == 'core'")
            .groupby("msp_name", as_index=False)
            .head(self.core_size)
        )
        return set(msp_content["gene_id"])

    def _streaming_counter(
        self, header: AlignmentHeader, ref_json: dict
    ) -> tuple[StreamingCounter, Path | None]:
        """Create a StreamingCounter (and the unsorted file of filtered alignments)

        The unsorted filtered alignments are written as uncompressed BAM: an
        unsorted CRAM would fetch reference sequences in random order.
        """
        strain_unsorted = None
        strain_genes = None
        if self.keep_filtered_alignments:
            strain_genes = self.get_strain_gene_ids(ref_json)
            strain_unsorted = Path(mkstemp(dir=self.meteor.tmp_dir, suffix=".bam")[1])
        streamer = StreamingCounter(
            header,
            self.counting_type,
            self.identity_threshold,
            strain_unsorted,
            strain_genes,
        )
        return streamer, strain_unsorted

    def _finish_streaming(
        self,
        streamer: StreamingCounter,
        strain_unsorted: Path | None,
        cramfile_strain: Path,
        count_file: Path,
        stage1_json_data: dict,
        stage1_json: Path,
        ref_json: dict,
    ) -> None:
        """Write the count table, the census and the sorted strain CRAM"""
        start = perf_counter()
        abundance = streamer.finish()
        if streamer.strain_out is not None:
            streamer.strain_out.close()
        logging.info("Counted reads in %f seconds", perf_counter() - start)
        start = perf_counter()
        self.write_stat(count_file, abundance, streamer.database)
        logging.info("Count table written in %f seconds", perf_counter() - start)
        total_read_count = stage1_json_data["mapping"]["total_read_count"]
        config = self.set_counter_config(
            total_read_count, streamer.counted_reads, count_file
        )
        stage1_json_data.update(config)
        self.save_config(stage1_json_data, stage1_json)
        if strain_unsorted is not None:
            start = perf_counter()
            self.sort_index_strain(strain_unsorted, cramfile_strain, ref_json)
            logging.info(
                "Filtered alignments sorted in %f seconds", perf_counter() - start
            )
        else:
            logging.info(
                "Cram file is not kept (--kf). Strain analysis will require a new mapping."
            )

    def sort_index_strain(
        self, unsorted: Path, cramfile_strain: Path, ref_json: dict
    ) -> None:
        reference = (
            self.meteor.ref_dir
            / ref_json["reference_file"]["fasta_dir"]
            / ref_json["reference_file"]["fasta_filename"]
        )
        sort(
            "-o",
            str(cramfile_strain.resolve()),
            "-@",
            str(self.meteor.threads),
            "-O",
            "cram",
            "--reference",
            str(reference.resolve()),
            str(unsorted.resolve()),
            catch_stdout=False,
        )
        index(str(cramfile_strain.resolve()))
        unsorted.unlink(missing_ok=True)

    def launch_mapping_and_counting(
        self,
        cramfile_strain: Path,
        count_file: Path,
        ref_json: dict,
        stage1_json: Path,
    ) -> None:
        """Map with bowtie2 and count in the same pass (no raw CRAM round trip).

        The raw CRAM is only written when --ka is requested.
        """
        state: dict = {}

        def factory(header: AlignmentHeader) -> StreamingCounter:
            streamer, strain_unsorted = self._streaming_counter(header, ref_json)
            state["streamer"] = streamer
            state["strain_unsorted"] = strain_unsorted
            return streamer

        pysam.set_verbosity(0)
        self._build_mapper(self.keep_all_alignments, factory).execute()
        stage1_json_data = self.read_json(stage1_json)
        start = perf_counter()
        self._finish_streaming(
            state["streamer"],
            state["strain_unsorted"],
            cramfile_strain,
            count_file,
            stage1_json_data,
            stage1_json,
            ref_json,
        )
        logging.info("Completed count table writing in %f seconds", perf_counter() - start)

    def set_counter_config(self, total_read_count: int, counted_reads: float, count_file: Path) -> dict:
        """Save in the json essential info
        :param total_read_count: (int) Number of reads mapped on the catalogue
        :param counted_read: (float) Number of reads counted after filtering
        :param count_file: (Path) Count file
        :return: (Dict) dictionnary data
        """
        return {
            "counting": {
                "counted_reads": counted_reads,
                "final_mapping_rate": round(100 * counted_reads / total_read_count, 2),
                "identity_threshold": round(self.identity_threshold, 2),
                "count_file": count_file.name,
            }
        }

    def write_stat(
        self, output: Path, abundance: dict[int, int | float], database: dict[int, int]
    ) -> None:
        """Write count table.

        :param output: [STRING] = output filename
        :param abundance_dict: [DICT] = contains abundance of each reference gene
        :param database: [DICT] = contains length of each reference gene
        """
        with lzma.open(output, "wt", preset=0) as out:
            out.write("gene_id\tgene_length\tvalue\n")
            for gene in sorted(abundance.keys()):
                out.write(f"{gene}\t{database[gene]}\t{abundance[gene]}\n")

    def launch_counting(
        self,
        raw_cramfile: Path,
        cramfile_strain: Path,
        count_file: Path,
        ref_json: dict,
        stage1_json_data: dict,
        stage1_json: Path,
    ):
        """Function that count reads from a cram file, using the given methods in count:
        "total" or "shared" or "unique".

        :param raw_cramfile: (Path) A path to the input cram file
        :param cramfile_strain: (Path) A path to the output cram file for strain analysis
        :param count_file: (Path) A path to the output count file
        """
        if not raw_cramfile.exists():
            logging.error(
                "Cram file %s is not available to perform a new counting. Please consider to re-map with --ka option.",
                raw_cramfile,
            )
            sys.exit(1)
        else:
            logging.info("Launch counting")
        pysam.set_verbosity(0)
        with AlignmentFile(
            str(raw_cramfile.resolve()), threads=self.meteor.threads
        ) as cramdesc:
            streamer, strain_unsorted = self._streaming_counter(
                cramdesc.header, ref_json
            )
            feed = streamer.feed
            for element in cramdesc:
                feed(element)
        self._finish_streaming(
            streamer,
            strain_unsorted,
            cramfile_strain,
            count_file,
            stage1_json_data,
            stage1_json,
            ref_json,
        )

    def execute(self) -> None:
        """Compute the mapping"""
        mapping_done = True
        ref_json = self.get_reference_info(self.meteor.ref_dir)
        Component.check_catalogue(ref_json)
        self.meteor.ref_name = ref_json["reference_info"]["reference_name"]
        if self.identity_user:
            self.identity_threshold = self.identity_user
        elif ref_json["reference_info"]["database_type"] == "complete":
            self.identity_threshold = self.DEFAULT_IDENTITY_THRESHOLD_COMPLETE
        else:
             self.identity_threshold = self.DEFAULT_IDENTITY_THRESHOLD_TAXO

        census_json_files = list(
            self.meteor.fastq_dir.glob("*_census_stage_0.json")
        )
        if len(census_json_files) == 0:
            logging.error(
                "No *_census_stage_0.json file found in %s",
                self.meteor.fastq_dir,
            )
            sys.exit(1)

        #  mapping of each sample against reference
        for library in census_json_files:
            census_json = self.read_json(library)
            sample_name = census_json["sample_info"]["sample_name"]
            stage1_dir = self.meteor.mapping_dir / sample_name
            stage1_dir.mkdir(exist_ok=True, parents=True)
            stage1_json = (
                self.meteor.mapping_dir
                / sample_name
                / f"{sample_name}_census_stage_1.json"
            )
            self.json_data[library] = {
                "census": census_json,
                "directory": stage1_dir,
                "Stage1FileName": stage1_json,
                "reference": ref_json,
            }
            if not stage1_json.exists():
                mapping_done = False
        # mapping already done and no overwriting
        cram_file = stage1_dir / f"{sample_name}.cram"
        count_file = stage1_dir / f"{sample_name}.tsv.xz"
        start = perf_counter()
        if mapping_done:
            logging.info(
                "Mapping already done for sample: %s",
                sample_name,
            )
            logging.info("Skipped !")
            # running counter from the raw cram
            stage1_json_data = self.read_json(stage1_json)
            raw_cram_file = stage1_dir / stage1_json_data["mapping"]["mapping_file"]
            self.launch_counting(
                raw_cram_file,
                cram_file,
                count_file,
                ref_json,
                stage1_json_data,
                stage1_json,
            )
            logging.info("Completed counting in %f seconds", perf_counter() - start)
        else:
            logging.info("Launch mapping and counting")
            self.launch_mapping_and_counting(cram_file, count_file, ref_json, stage1_json)
            logging.info(
                "Completed mapping and counting in %f seconds", perf_counter() - start
            )
            raw_cram_file = (
                stage1_dir / self.read_json(stage1_json)["mapping"]["mapping_file"]
            )
        if not self.keep_all_alignments:
            logging.info(
                "Raw cram file is not kept (--ka). "
                "Re-counting operation will need to be performed from scratch."
            )
            raw_cram_file.unlink(missing_ok=True)
            raw_cram_file.with_suffix(".cram.crai").unlink(missing_ok=True)

        # Not sure if it's a good idea to delete a temporary file
        rmtree(self.meteor.tmp_dir, ignore_errors=True)
