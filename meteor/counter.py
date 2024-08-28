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
import pandas as pd
from dataclasses import dataclass, field
from tempfile import mkdtemp, mkstemp
from pathlib import Path
from meteor.mapper import Mapper
from meteor.session import Session, Component
from typing import Generator, Iterator, ClassVar
from collections import defaultdict
from itertools import chain
from pysam import index, idxstats, AlignmentFile, sort, AlignedSegment  # type: ignore[attr-defined]
from time import perf_counter
from shutil import rmtree


@dataclass
class Counter(Session):
    """Counter session map and count"""

    COUNTING_TYPES: ClassVar[list[str]] = ["total", "smart_shared", "unique"]
    DEFAULT_COUNTING_TYPE: ClassVar[str] = "smart_shared"
    NO_IDENTITY_THRESHOLD: ClassVar[float] = 0.0
    DEFAULT_IDENTITY_THRESHOLD_COMPLETE: ClassVar[float] = 0.95
    DEFAULT_IDENTITY_THRESHOLD_TAXO: ClassVar[float] = 0.97

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
        list_fastq_path = []
        # loop on each library
        for (
            _,
            dict_data,
        ) in self.json_data.items():  # pylint: disable=unused-variable
            census = dict_data["census"]
            sample_file = census["sample_file"]
            # reindexing this library reads and fill FLibraryIndexerReport
            list_fastq_path += [str(self.meteor.fastq_dir / sample_file["fastq_file"])]
        else:
            # mapping this library on the reference
            mapping_process = Mapper(
                self.meteor,
                dict_data,
                list_fastq_path,
                self.mapping_type,
                self.trim,
                self.alignment_number,
                self.identity_threshold,
            )
            mapping_process.execute()

    def write_table(self, cramfile: Path, outfile: Path) -> None:
        """Function that create a count table using pysam. First index the cram file,
        then count reads using the function idxstats from pysam, and output a count
        table.

        :param cramfile: (Path) cram file to count
        :param outfile: (Path) count table
        """
        # indStats = pd.read_csv(StringIO(pysam.idxstats(cramfile)), sep = '\t',
        # header = None, names = ['contig', 'length', 'mapped', 'unmapped'])
        # create count table
        table = idxstats(str(cramfile.resolve()))
        # write the count table
        with lzma.open(outfile, "wt", preset=0) as out:
            out.write("gene_id\tgene_length\tvalue\n")
            for line in table.split("\n")[:-2]:
                # print(int(line.split("\t")[0]))
                # print(unique_on_gene[int(line.split("\t")[0])])
                # print(line.split("\t")[2])
                # assert unique_on_gene[int(line.split("\t")[0])] == int(
                #     line.split("\t")[2]
                # ), "Gene length does not match counted nucleotides"
                s = "\t".join(line.split("\t")[0:3])
                out.write(f"{s}\n")

    def get_aligned_nucleotides(self, element) -> Iterator[int]:
        """Select aligned nucleotides
        :param element: Alignment object
        :return: Aligned item
        """
        yield from (item[1] for item in element.cigartuples if item[0] < 3)

    def filter_alignments(
        self, cramdesc: AlignmentFile
    ) -> tuple[dict[str, list[AlignedSegment]], dict[str, list[int]]]:
        """Filter read according to their identity with reference and reads with multiple
        alignments with different score. We keep the best scoring reads when total count is
        applied.
        Shared count keeps multiple alignments when score are equal

        :param cramdesc [STR] = CRAM file to count
        :return: A tuple with database [DICT] = contains length of reference genes.
                                        key :
                                        value : reference gene
                And  genes [DICT] = contains the set of reference genes.
                                    key : read_id
        """
        tmp_score: dict[str, float] = {}
        genes: dict[str, list[int]] = {}
        # contains a list of alignment of each read
        reads: dict[str, list[AlignedSegment]] = {}
        for element in cramdesc:
            assert element.query_name is not None and element.reference_name is not None

            # identity = (element.query_length - element.get_tag("NM")) / element.query_length
            # identity = 1.0 - (element.get_tag("NM") / element.query_alignment_length)
            ali = sum(self.get_aligned_nucleotides(element))
            # identity = (element.query_length - element.get_tag("NM")) / element.query_length
            identity = (ali - int(element.get_tag("NM"))) / ali
            # if lower than the identity threshold
            # we ignore the read
            if identity < self.identity_threshold:
                continue
            # Only if we use score
            # if not element.has_tag("AS"):
            #     raise ValueError("Missing 'AS' field.")
            read_id: str = element.query_name
            # print(read_id, element.query_alignment_length)
            # get alignment score
            # Meteor do not take in account the alignement score
            # But the identity
            score = identity
            # score = element.get_tag("AS")
            # get previous score for the read
            prev_score = tmp_score.get(read_id, 0)
            # higher = more similar
            # https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
            # With 1.0 identity, score is 0
            # With 0.975 identity, score is -16...
            # if same score
            if prev_score == score:
                # add the genes to the list if it doesn't exist
                reads[read_id].append(element)
                genes[read_id].append(int(element.reference_name))
            # case new score is higher
            elif prev_score < score:
                # set the new score
                tmp_score[read_id] = score
                # We keep the new score and forget the previous one
                reads[read_id] = [element]
                genes[read_id] = [int(element.reference_name)]
        return reads, genes

    def uniq_from_mult(
        self, reads: dict, genes: dict, database: dict
    ) -> tuple[defaultdict, dict, dict]:
        """
        Function that filter unique reads from all reads. Multiple reads are
        reads that map to more than one genes. And Unique reads are reads that map
        only on one genes.

        :params reads: [DICT] = dictionary containing reads as key and pysam alignment as value
        :params genes: [DICT] = dictionary containing reads as key and a list of
                                genes where read mapped as value
        :params database: [DICT] = contains for each reference genes the number of
                                unique reads
        :return: A tuple with:
            genes [DICT] = the dictionary without unique reads
            unique [DICT] = nb of unique read of each reference genes
        """
        unique_reads: defaultdict[
            str, list[pysam.libcalignedsegment.AlignedSegment]
        ] = defaultdict(list)
        unique_reads_list = []
        unique_on_gene: dict[int, int] = dict.fromkeys(database, 0)
        for read_id in genes:
            # if read map to several gene
            if len(genes[read_id]) != 1:
                continue
            # otherwise read map with best score against only one genes
            # we keep read
            unique_reads[read_id] = reads[read_id]
            # get genes id
            gene = genes[read_id][0]
            # add to unique read dictionnary
            unique_on_gene[gene] += 1
            unique_reads_list.append(read_id)
        # delete unique reads
        for read_id in unique_reads_list:
            # delete the read_id from TMP dictionnary
            del genes[read_id]
        return unique_reads, genes, unique_on_gene

    def compute_co(self, genes_mult: dict, unique_on_gene: dict) -> tuple[dict, dict]:
        """Compute genes specific coefficient "Co" for each multiple read.

        :param genes_mult: [DICT] = Contains for each multiple read, all the reference
                                genes name where it mapped
        :param unique_on_gene: [DICT] = nb of unique read of each reference genes
        :return: Co [DICT] = contains coefficient values for each couple (read, genes)
                read_dict [DICT] = contains all the multiple reads of a genes
                                read_id : genes
                                value : list a multiple reads
        """
        co_dict: dict[tuple[str, int], float] = {}
        read_dict: dict[str, list[str]] = {}
        for read_id in genes_mult:
            # for a multiple read, gets nb of unique reads from each genes
            # total number of unique reads
            som = sum(
                unique_on_gene[genes] for genes in genes_mult[read_id]
            )  # Sum the values directly
            # som = reduce(lambda x, y: x + y, s)
            # No unique counts on these genes
            # If no unique counts:
            if som == 0:
                # Specific count of Meteor
                # 1 / nb genes aligned by the read
                for genes in genes_mult[read_id]:
                    co_dict[(read_id, genes)] = 1.0 / len(genes_mult[read_id])
                    read_dict.setdefault(genes, []).append(read_id)
                # Normally we continue here
                continue
            # We get the unique set of duplicated genes
            # These genes are mapped several time
            # duplicated_genes = set([genes for genes in genes_mult[read_id] if genes_mult[read_id].count(genes) > 1])
            duplicated_genes = {
                genes
                for genes in genes_mult[read_id]
                if genes_mult[read_id].count(genes) > 1
            }
            # otherwise
            for genes in genes_mult[read_id]:
                # get the nb of unique reads
                nb_unique = unique_on_gene[genes]
                if nb_unique == 0:
                    # test if meteor do that
                    # co_dict[(read_id, genes)] = 1.0 / len(genes_mult[read_id])
                    continue
                # else:
                # calculate Co of multiple read for the given genes
                if genes in duplicated_genes:
                    if (read_id, genes) in co_dict:
                        co_dict[(read_id, genes)] += nb_unique / float(som)
                    else:
                        co_dict[(read_id, genes)] = nb_unique / float(som)
                else:
                    co_dict[(read_id, genes)] = nb_unique / float(som)
                # append to the dict
                read_dict.setdefault(genes, []).append(read_id)
        # get all the multiple reads of a genes and uniquify
        for genes, reads in read_dict.items():
            read_dict[genes] = list(set(reads))
        return read_dict, co_dict

    def get_co_coefficient(self, gene, read_dict: dict, co_dict: dict) -> Generator:
        """
        Compute genes specific coefficient "Co" for each multiple read.

        :param gene: [str] = A gene name
        :param read_dict: [DICT] = contains all the multiple reads of a genes
        :param co_dict: [DICT] = contains coefficient values for each couple (read, genes)
        :return: [float] = coefficient obtain for each pair (read, gene)
        """
        for read in read_dict[gene]:
            yield co_dict[(read, gene)]

    def compute_abm(self, read_dict: dict, co_dict: dict, database: dict) -> dict:
        """Compute multiple reads abundance for each genes.
        Multiple reads abundance of a genes is equal to the sum of all Co
        coefficient

        :param read_dict: [DICT] = list of multiple read for each reference genes
        :param co_dict: [DICT] = contains coefficient values for each couple (read, genes)
        :param multiple_dict: [DICT] = abundance of multiple reads for each genes
        :return: multiple_dict [DICT] = abundance of multiple reads for each genes
        """
        multiple_dict = dict.fromkeys(database, 0)
        for gene in read_dict:
            # get a list of all the Co coefficient for each reference genes
            # sum them, (we obtain the abundance)
            multiple_dict[gene] = sum(self.get_co_coefficient(gene, read_dict, co_dict))
        return multiple_dict

    def compute_abs_meteor(
        self, database: dict, unique_dict: dict, multiple_dict: dict
    ) -> dict:
        """Compute the abundance of a each genes.
            Abundance = Abundance unique reads + Abundance multiple reads

        :param unique_dict: [DICT] = nb of unique read of each reference genes
        :param multiple_dict: [DICT] = abundance of multiple reads for each genes
        :return: abundance_dict [DICT] = contains abundance of each reference genes
        """
        return {gene: unique_dict[gene] + multiple_dict[gene] for gene in database}

    def compute_abs(self, unique_dict: dict, multiple_dict: dict) -> dict:
        """Compute the abundance of a each genes.
            Abundance = Abundance unique reads + Abundance multiple reads

        :param unique_dict: [DICT] = nb of unique read of each reference genes
        :param multiple_dict: [DICT] = abundance of multiple reads for each genes
        :return: abundance_dict [DICT] = contains abundance of each reference genes
        """
        return {
            genes: unique_dict[genes] + multiple_dict[genes] for genes in unique_dict
        }

    def write_stat(self, output: Path, abundance_dict: dict, database: dict) -> None:
        """Write count table.

        :param output: [STRING] = output filename
        :param abundance_dict: [DICT] = contains abundance of each reference genes
        :param database: [DICT] = contrains length of each reference genes
        """
        with lzma.open(output, "wt", preset=0) as out:
            out.write("gene_id\tgene_length\tvalue\n")
            for genes, abundance in sorted(abundance_dict.items()):
                out.write(f"{genes}\t{database[genes]}\t{abundance}\n")

    def save_cram_strain(
        self,
        outcramfile: Path,
        cramdesc: AlignmentFile,
        read_chain: chain,
        ref_json: dict,
    ):
        """Writing the filtered CRAM file for strain analysis.

        :param outsamfile: [Path] Temporary cram file
        :param cramdesc: Pysam cram descriptor
        :param read_chain: [chain] Chain of pysam reads objects
        """
        reference = (
            self.meteor.ref_dir
            / ref_json["reference_file"]["fasta_dir"]
            / ref_json["reference_file"]["fasta_filename"]
        )
        msp_file = (
            self.meteor.ref_dir
            / ref_json["reference_file"]["database_dir"]
            / ref_json["annotation"]["msp"]["filename"]
        )
        msp_content = (
            pd.read_csv(
                msp_file,
                sep="\t",
                # names=["msp_name", "gene_id", "gene_name", "gene_category"],
                header=0,
            )
            .query("gene_category == 'core'")
            .groupby("msp_name", as_index=False)
            .head(self.core_size)
        )
        # Sort the 'gene_id' values and convert directly to a list
        gene_id_set = set(msp_content["gene_id"])
        with AlignmentFile(
            str(outcramfile.resolve()),
            "wc",
            template=cramdesc,
            reference_filename=str(reference.resolve()),
            threads=self.meteor.threads,
        ) as total_reads:
            for element in read_chain:
                if int(element.reference_name) in gene_id_set:
                    # if int(element.reference_name) in ref_json["reference_file"]:
                    total_reads.write(element)

    def launch_counting(
        self,
        raw_cramfile: Path,
        cramfile_strain: Path,
        count_file: Path,
        ref_json: dict,
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
            # create a dictionary containing the length of reference genes
            # get name of reference sequence
            references = [int(ref) for ref in cramdesc.references]
            # get reference length
            lengths = cramdesc.lengths
            # All references with their length
            database = dict(zip(references, lengths))
            # Filter reads according to quality and score
            reads, genes = self.filter_alignments(cramdesc)
            # Filter reads according to quality, score and counting
            # genes mapped by multiple reads
            unique_reads, genes_mult, unique_on_gene = self.uniq_from_mult(
                reads, genes, database
            )
            if self.counting_type == "smart_shared":
                # For multiple reads compute Co
                read_dict, coef_read = self.compute_co(genes_mult, unique_on_gene)
                # calculate abundance
                multiple = self.compute_abm(read_dict, coef_read, database)
                # Calculate reference abundance & write count table
                abundance = self.compute_abs_meteor(database, unique_on_gene, multiple)
                # abundance = self.compute_abs(unique_on_gene, multiple)
                self.write_stat(count_file, abundance, database)

            else:
                if self.counting_type == "unique":
                    reads = unique_reads
                cramfile_unsorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                cramfile_sorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                # self.json_data[library]["directory"]
                #     / f"{sample_info['sample_name']}_raw.cram"
                # self.save_cram(
                #     cramfile_unsorted, cramdesc, list(chain(*reads.values())), ref_json
                # )
                reference = (
                    self.meteor.ref_dir
                    / ref_json["reference_file"]["fasta_dir"]
                    / ref_json["reference_file"]["fasta_filename"]
                )
                # Save a temporary cram for the counting
                with AlignmentFile(
                    str(cramfile_unsorted.resolve()),
                    "wc",
                    template=cramdesc,
                    reference_filename=str(reference.resolve()),
                    threads=self.meteor.threads,
                ) as total_reads:
                    for element in chain.from_iterable(reads.values()):
                        # if int(element.reference_name) in ref_json["reference_file"]:
                        total_reads.write(element)
                sort(
                    "-o",
                    str(cramfile_sorted.resolve()),
                    "-@",
                    str(self.meteor.threads),
                    "-O",
                    "cram",
                    str(cramfile_unsorted.resolve()),
                    catch_stdout=False,
                )
                self.write_table(cramfile_sorted, count_file)
            if self.keep_filtered_alignments:
                cramfile_strain_unsorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                self.save_cram_strain(
                    cramfile_strain_unsorted,
                    cramdesc,
                    chain.from_iterable(reads.values()),
                    ref_json,
                )
                sort(
                    "-o",
                    str(cramfile_strain.resolve()),
                    "-@",
                    str(self.meteor.threads),
                    "-O",
                    "cram",
                    str(cramfile_strain_unsorted.resolve()),
                    catch_stdout=False,
                )
                index(str(cramfile_strain.resolve()))
            else:
                logging.info(
                    "Cram file is not kept (--kf). Strain analysis will require a new mapping."
                )

    def execute(self) -> None:
        """Compute the mapping"""
        mapping_done = True
        try:
            # Get the ini ref
            ref_json_file_list = list(self.meteor.ref_dir.glob("**/*_reference.json"))
            assert len(ref_json_file_list) == 1
            ref_json_file = ref_json_file_list[0]
            ref_json = self.read_json(ref_json_file)
            self.meteor.ref_name = ref_json["reference_info"]["reference_name"]
            if not self.identity_user:
                if ref_json["reference_info"]["database_type"] == "complete":
                    self.identity_threshold = self.DEFAULT_IDENTITY_THRESHOLD_COMPLETE
                else:
                    self.identity_threshold = self.DEFAULT_IDENTITY_THRESHOLD_TAXO
            else:
                self.identity_threshold = self.identity_user
        except AssertionError:
            logging.error(
                "No *_reference.json file found in %s. "
                "One *_reference.json is expected",
                self.meteor.ref_dir,
            )
            sys.exit(1)
        try:
            census_json_files = list(
                self.meteor.fastq_dir.glob("*_census_stage_0.json")
            )
            assert len(census_json_files) > 0

            #  mapping of each sample against reference
            for library in census_json_files:
                census_json = self.read_json(library)
                sample_info = census_json["sample_info"]
                stage1_dir = self.meteor.mapping_dir / sample_info["sample_name"]
                stage1_dir.mkdir(exist_ok=True, parents=True)
                # if self.pysam_test:
                self.json_data[library] = {
                    "census": census_json,
                    "directory": stage1_dir,
                    "Stage1FileName": stage1_dir
                    / f"{sample_info['sample_name']}_census_stage_1.json",
                    "reference": ref_json,
                }
                if not self.json_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            else:
                # mapping already done and no overwriting
                if mapping_done:
                    logging.info(
                        "Mapping already done for sample: %s",
                        sample_info["sample_name"],
                    )
                    logging.info("Skipped !")
                else:
                    logging.info("Launch mapping")
                    self.launch_mapping()
                # running counter
                raw_cram_file = (
                    self.json_data[library]["directory"]
                    / f"{sample_info['sample_name']}_raw.cram"
                )
                cram_file = (
                    self.json_data[library]["directory"]
                    / f"{sample_info['sample_name']}.cram"
                )
                count_file = (
                    self.json_data[library]["directory"]
                    / f"{sample_info['sample_name']}.tsv.xz"
                )
                start = perf_counter()
                self.launch_counting(raw_cram_file, cram_file, count_file, ref_json)
                logging.info("Completed counting in %f seconds", perf_counter() - start)
                if not self.keep_all_alignments:
                    logging.info(
                        "Raw cram file is not kept (--ka). "
                        "Re-counting operation will need to be performed from scratch."
                    )
                    raw_cram_file.unlink(missing_ok=True)
                    raw_cram_file.with_suffix(".cram.crai").unlink(missing_ok=True)
        except AssertionError:
            logging.error(
                "No *_census_stage_0.json file found in %s",
                self.meteor.fastq_dir,
            )
            sys.exit(1)
        else:
            # Not sure if it's a good idea to delete a temporary file
            rmtree(self.meteor.tmp_dir, ignore_errors=True)
