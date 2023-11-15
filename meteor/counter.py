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

from functools import reduce
import gzip
import bz2
import lzma
import logging
import sys
import pysam
from dataclasses import dataclass, field
from tempfile import mkdtemp, mkstemp, NamedTemporaryFile
import tempfile
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
from meteor.mapper import Mapper
from meteor.session import Session, Component
from typing import Type, Dict, Generator, List, Tuple, Iterator
from collections import defaultdict
from itertools import chain
from pysam import index, idxstats, AlignmentFile, sort, AlignedSegment  # type: ignore[attr-defined]
from time import perf_counter
from shutil import rmtree, copy


@dataclass
class Counter(Session):
    """Counter session map and count"""

    meteor: Type[Component]
    counting_type: str
    mapping_type: str
    trim: int
    identity_threshold: float
    alignment_number: int
    counting_only: bool
    mapping_only: bool
    keep_sam: bool = False
    keep_bam: bool = False
    # pysam_test: bool = True
    ini_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.meteor.tmp_path:
            self.meteor.tmp_path.mkdir(exist_ok=True)
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.mapping_dir.mkdir(exist_ok=True)

    def count_index_fastq(
        self, fastq_file: Path, output_desc: tempfile._TemporaryFileWrapper
    ) -> tuple:
        """Count the number of bases

        :param fastq_file: Path object of the fastq_file
        :param output_desc: Output tempfile._TemporaryFileWrapper descriptor
        :return: (Tuple) A tuple giving the read count and base count
        """
        read_count = 0
        base_count = 0
        if fastq_file.suffix == ".gz":
            in_fq = gzip.open(fastq_file, "rt")
        elif fastq_file.suffix == ".bz2":
            in_fq = bz2.open(fastq_file, "rt")
        elif fastq_file.suffix == ".xz":
            in_fq = lzma.open(fastq_file, "rt")
        else:
            in_fq = open(fastq_file, "rt", encoding="UTF-8")
        # read input fastq line by line
        with in_fq:
            for read_count, line in enumerate(
                in_fq, start=1
            ):  # pylint: disable=unused-variable
                output_desc.write(f"@{read_count}\n")
                # read the sequence
                read = next(in_fq)
                output_desc.write(f"{read}")
                base_count += len(read.strip())
                # pass the plus
                output_desc.write(next(in_fq))
                # pass the quality
                output_desc.write(next(in_fq))
            output_desc.flush()
        return read_count, base_count

    def set_workflow_config(self, ref_ini) -> ConfigParser:
        """Write configuration file for reference genes

        :return: (ConfigParser) A configparser object
        """
        config = ConfigParser()
        config["worksession"] = {
            "meteor.reference.dir": self.meteor.ref_dir.name,
            "meteor.mapping.program": "bowtie2",
            "meteor.mapping.file.format": "bam",
            "meteor.cpu.count": str(self.meteor.threads),
        }
        config["main_reference"] = {
            "meteor.reference.name": self.meteor.ref_name,
            "meteor.matches": str(self.alignment_number),
            "meteor.mismatches": "5",
            "meteor.is.perc.mismatches": "1",
            "meteor.bestalignment": "1",
            "meteor.mapping.prefix.name": f"mapping_vs_{ref_ini['reference_info']['reference_name']}",
            "meteor.counting.prefix.name": f"vs_{ref_ini['reference_info']['reference_name']}",
        }
        return config

    def launch_mapping(self) -> None:
        """Create temporary indexed files and map against"""
        logging.info("Launch mapping")
        list_fastq_path = []
        # loop on each library
        for (
            library,
            dict_data,
        ) in self.ini_data.items():  # pylint: disable=unused-variable
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
                self.counting_type,
            )
            if not mapping_process.execute():
                raise ValueError("Error, TaskMainMapping failed")

    def write_table(self, bamfile: Path, outfile: Path) -> bool:
        """Function that create a count table using pysam. First index the BAM file,
        then count reads using the function idxstats from pysam, and output a count
        table.

        :param bamfile: (Path) BAM file to count
        :param outfile: (Path) count table
        """
        if not bamfile.with_suffix(".bam.bai").exists():
            # index the bam file
            index(str(bamfile.resolve()))
        # indStats = pd.read_csv(StringIO(pysam.idxstats(bamfile)), sep = '\t',
        # header = None, names = ['contig', 'length', 'mapped', 'unmapped'])
        # create count table
        table = idxstats(str(bamfile.resolve()))
        # write the count table
        with outfile.open("wt", encoding="UTF-8") as out:
            out.write(
                f"{self.meteor.gene_column}\t{self.meteor.gene_length_column}\t{self.meteor.value_column}\n"
            )
            for line in table.split("\n")[:-2]:
                s = "\t".join(line.split("\t")[0:3])
                out.write(f"{s}\n")
        return True

    def get_aligned_nucleotides(self, element) -> Iterator[int]:
        """Select aligned nucleotides
        :param element: Alignment object
        :return: Aligned item
        """
        for item in element.cigartuples:
            # I don't know this cigar yet
            if item[0] < 3:
                yield item[1]

    def filter_alignments(
        self, samdesc: AlignmentFile
    ) -> tuple[defaultdict, defaultdict]:
        """Filter read according to their identity with reference and reads with multiple
        alignments with different score. We keep the best scoring reads when total count is
        applied.
        Shared count keeps multiple alignments when score are equal

        :param samdesc [STR] = SAM file to count
        :return: A tuple with database [DICT] = contains length of reference genes.
                                        key :
                                        value : reference gene
                And  genes [DICT] = contains the set of reference genes.
                                    key : read_id
        """
        tmp_score: Dict[str, float] = {}
        genes: defaultdict[str, List[int]] = defaultdict(list)
        # contains a list of alignment of each read
        reads: defaultdict[str, List[AlignedSegment]] = defaultdict(list)
        for element in samdesc:
            # identity = (element.query_length - element.get_tag("NM")) / element.query_length
            # identity = 1.0 - (element.get_tag("NM") / element.query_alignment_length)
            ali = sum(list(self.get_aligned_nucleotides(element)))
            # identity = (element.query_length - element.get_tag("NM")) / element.query_length
            identity = (ali - int(element.get_tag("NM"))) / ali
            # if lower than the identity threshold
            # we ignore the read
            if identity < self.identity_threshold:
                continue
            # Only if we use score
            # if not element.has_tag("AS"):
            #     raise ValueError("Missing 'AS' field.")
            read_id: str | None = element.query_name
            if not read_id:
                continue
            # print(read_id, element.query_alignment_length)
            # get alignment score
            # Meteor do not take in account the alignement score
            # But the identity
            score = identity
            # score = element.get_tag("AS")
            # get previous score for the read
            prev_score = tmp_score.get(read_id, score)
            # higher = more similar
            # https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#scores-higher-more-similar
            # With 1.0 identity, score is 0
            # With 0.975 identity, score is -16...
            # if same score
            if prev_score == score:
                tmp_score[read_id] = score
                # add the genes to the list if it doesn't exist
                # Meteor uncomment if we don't want two counts for the same gene
                # if int(element.reference_name) not in set(genes[read_id]):
                # if int(element.reference_name) == 182 and int(element.reference_name) in set(genes[read_id]):
                #     print("gene mapped several time")
                if not element.reference_name:
                    continue
                reads[read_id].append(element)
                genes[read_id].append(int(element.reference_name))
                # else:
                #     print(reads[read_id])
            # case new score is higher
            elif prev_score < score:
                # set the new score
                tmp_score[read_id] = score
                # In best counting, it happens because several alignment can be selected by bowtie
                # if self.counting_type == "best":
                #     raise ValueError("Bam file contains read aligned against multiple genes. "
                #                      "It should not, we aligned with the bowtie2 -k 1 option.")
                # In smart_shared and unique case, we reinitialise all
                # We keep the new score and forget the previous one
                # elif self.counting_type in ("unique", "smart_shared"):
                if not element.reference_name:
                    continue
                print(element)
                reads[read_id] = [element]
                genes[read_id] = [int(element.reference_name)]
                # In total count, we keep all but not meteor which performs an ex-aequo count
                # Case total
                # else:
                #     reads[read_id].append(element)
                #     genes[read_id].append(element.reference_name)
        # we need to check this following part
        # May be useless
        # if self.counting_type == "smart_shared":
        #     # We keep a unique reference per read
        #     # genes dict is only used for smart shared counting
        #     for read_id, reference_names in genes.items():
        #         genes[read_id] = list(set(reference_names))
        return reads, genes

    def uniq_from_mult(
        self, reads: defaultdict, genes: defaultdict, database: dict
    ) -> tuple[defaultdict, defaultdict, dict]:
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
            str, List[pysam.libcalignedsegment.AlignedSegment]
        ] = defaultdict(list)
        unique_reads_list = []
        # dict[str]:int
        # gene : 0
        # count = 0
        unique_on_gene: Dict[int, int] = dict.fromkeys(database, 0)
        # print("ici")
        # print(genes['1791'])
        for read_id in genes:
            # gene = genes[read_id][0]
            # if read map to several gene
            if len(genes[read_id]) != 1:
                continue
            # print(genes[read_id])
            # otherwise read map with best score against only one genes
            # we keep read
            unique_reads[read_id] = reads[read_id]
            # get genes id
            gene = genes[read_id][0]
            # if gene == 11208:
            #     print(read_id)
            #     count +=1
            #     print(f"I have unique {count}")
            # add to unique read dictionnary
            unique_on_gene[gene] += 1
            unique_reads_list.append(read_id)
        # delete unique reads
        for read_id in unique_reads_list:
            # delete the read_id from TMP dictionnary
            del genes[read_id]
        return unique_reads, genes, unique_on_gene

    def compute_co(
        self, genes_mult: defaultdict, unique_on_gene: dict
    ) -> tuple[Dict, Dict]:
        """Compute genes specific coefficient "Co" for each multiple read.

        :param genes_mult: [DICT] = Contains for each multiple read, all the reference
                                genes name where it mapped
        :param unique_on_gene: [DICT] = nb of unique read of each reference genes
        :return: Co [DICT] = contains coefficient values for each couple (read, genes)
                read_dict [DICT] = contains all the multiple reads of a genes
                                read_id : genes
                                value : list a multiple reads
        """
        co_dict: Dict[Tuple[str, int], float] = {}
        read_dict: Dict[str, List[str]] = {}
        for read_id in genes_mult:
            # for a multiple read, gets nb of unique reads from each genes
            s = [unique_on_gene[genes] for genes in genes_mult[read_id]]
            # total number of unique reads
            som = reduce(lambda x, y: x + y, s)
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

    def write_stat(self, output: Path, abundance_dict: dict, database: dict):
        """Write count table.

        :param output: [STRING] = output filename
        :param abundance_dict: [DICT] = contains abundance of each reference genes
        :param database: [DICT] = contrains length of each reference genes
        """
        with open(output, "wt", encoding="UTF-8") as out:
            out.write(
                f"{self.meteor.gene_column}\t{self.meteor.gene_length_column}\t{self.meteor.value_column}\n"
            )
            for genes, abundance in sorted(abundance_dict.items()):
                out.write(f"{genes}\t{database[genes]}\t{abundance}\n")
        return True

    def save_bam(self, outbamfile: Path, samdesc: AlignmentFile, read_list: list):
        """Writing the filtered SAM file.

        :param outsamfile: [Path] Temporary bam file
        :param samdesc: Pysam sam descriptor
        :param read_list: [List] List of pysam reads objects
        """
        with AlignmentFile(
            str(outbamfile.resolve()), "wb", template=samdesc
        ) as total_reads:
            for element in read_list:
                total_reads.write(element)

    # def launch_counting(self, bam_file: Path, count_file: Path) -> bool:
    def launch_counting(self, sam_file: Path, count_file: Path) -> bool:
        """Function that count reads from a BAM file, using the given methods in count:
        "total" or "shared" or "unique".

        :param sam_file: (Path) A path to the input bam file
        :param count_file: (Path) A path to the output count file
        """
        logging.info("Launch counting")
        # "-t", str(self.meteor.tmp_dir) + "/"
        # open the BAM file
        # with AlignmentFile(str(bam_file.resolve()), "rb") as bamdesc:
        with AlignmentFile(str(sam_file.resolve())) as samdesc:
            # create a dictionary containing the length of reference genes
            # get name of reference sequence
            references = [int(ref) for ref in samdesc.references]
            # get reference length
            lengths = samdesc.lengths
            # All references with their length
            database = dict(zip(references, lengths))
            # Filter reads according to quality and score
            reads, genes = self.filter_alignments(samdesc)
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
                if self.keep_bam:
                    bamfile = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                    bamfile_sorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                    self.save_bam(bamfile, samdesc, list(chain(*reads.values())))
                    sort(
                        "-o",
                        str(bamfile_sorted.resolve()),
                        "-@",
                        str(self.meteor.threads),
                        "-O",
                        "bam",
                        str(bamfile.resolve()),
                        catch_stdout=False,
                    )
                    bamfile_sorted = Path(
                        copy(
                            str(bamfile_sorted.resolve()),
                            str(sam_file.resolve().with_suffix(".bam")),
                        )
                    )
                    index(str(bamfile_sorted.resolve()))
                return self.write_stat(count_file, abundance, database)
            else:
                if self.counting_type == "unique":
                    reads = unique_reads
                bamfile = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                bamfile_sorted = Path(mkstemp(dir=self.meteor.tmp_dir)[1])
                self.save_bam(bamfile, samdesc, list(chain(*reads.values())))
                sort(
                    "-o",
                    str(bamfile_sorted.resolve()),
                    "-@",
                    str(self.meteor.threads),
                    "-O",
                    "bam",
                    str(bamfile.resolve()),
                    catch_stdout=False,
                )
                if self.keep_bam:
                    bamfile_sorted = Path(
                        copy(
                            str(bamfile_sorted.resolve()),
                            str(sam_file.resolve().with_suffix(".bam")),
                        )
                    )
                    index(str(bamfile_sorted.resolve()))
                return self.write_table(bamfile_sorted, count_file)

    def execute(self) -> bool:
        """Compute the mapping"""
        mapping_done = True
        try:
            # Get the ini ref
            ref_ini_file_list = list(self.meteor.ref_dir.glob("**/*_reference.ini"))
            assert len(ref_ini_file_list) == 1
            ref_ini_file = ref_ini_file_list[0]
            ref_ini = ConfigParser()
            with open(ref_ini_file, "rt", encoding="UTF-8") as ref:
                ref_ini.read_file(ref)
            self.meteor.ref_name = ref_ini["reference_info"]["reference_name"]
        except AssertionError:
            logging.error(
                "Error, no *_reference.ini file found in %s. "
                "One *_reference.ini is expected",
                self.meteor.ref_dir,
            )
            sys.exit()
        try:
            census_ini_files = list(self.meteor.fastq_dir.glob("*_census_stage_0.ini"))
            assert len(census_ini_files) > 0
            #  MAPPING THIS LIBRARY ON MAIN REFERENCE
            for library in census_ini_files:
                census_ini = ConfigParser()
                with open(library, "rt", encoding="UTF-8") as cens:
                    census_ini.read_file(cens)
                    sample_info = census_ini["sample_info"]
                    # sample_file = census_ini['sample_file'] # reference
                    #
                    # if self.pysam_test:
                    stage1_dir = self.meteor.mapping_dir / sample_info["sample_name"]
                    # else:
                    # stage1_dir = (
                    #     self.meteor.mapping_dir / sample_info["sample_name"] /
                    #     f"mapping_vs_{self.meteor.ref_name}_{census_ini['sample_info']['full_sample_name']}"
                    # )
                    # sample_info["full_sample_name"]
                    stage1_dir.mkdir(exist_ok=True, parents=True)
                    # if self.pysam_test:
                    self.ini_data[library] = {
                        "census": census_ini,
                        "directory": stage1_dir,
                        "Stage1FileName": stage1_dir
                        / f"{sample_info['sample_name']}_census_stage_1.ini",
                        "reference": ref_ini,
                    }
                    # else:
                    # self.ini_data[library] = {
                    #     "census": census_ini,
                    #     "directory": stage1_dir,
                    #     "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1"),
                    #     "reference": ref_ini
                    # }
                if not self.ini_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            else:
                if not self.counting_only:
                    # mapping already done and no overwriting
                    if mapping_done:
                        logging.info(
                            "Mapping already done for sample: %s",
                            sample_info["sample_name"],
                        )
                        logging.info("Skipped !")
                    else:
                        self.launch_mapping()
                if not self.mapping_only:
                    logging.info("Launch mapping")
                    config = self.set_workflow_config(ref_ini)
                    workflow_ini = self.meteor.mapping_dir / "workflow.ini"
                    self.save_config(config, workflow_ini)
                    # bam_file = self.ini_data[library]["directory"] / f"{sample_info['sample_name']}.bam"
                    # test
                    sam_file = (
                        self.ini_data[library]["directory"]
                        / f"{sample_info['sample_name']}.sam"
                    )
                    count_file = (
                        self.ini_data[library]["directory"]
                        / f"{sample_info['sample_name']}.tsv"
                    )
                    # if self.pysam_test:
                    start = perf_counter()
                    # self.launch_counting(bam_file, count_file)
                    self.launch_counting(sam_file, count_file)
                    logging.info(
                        "Completed counting in %f seconds", perf_counter() - start
                    )
                    # else:
                    # self.launch_counting(workflow_ini)
                    if not self.keep_sam:
                        sam_file.unlink(missing_ok=True)
                        # bam_file.with_suffix(".bam.bai").unlink(missing_ok=True)
                        # self.ini_data[library]["Stage1FileName"].unlink(missing_ok=True)
            logging.info("Done ! Job finished without errors ...")
        except AssertionError:
            logging.error(
                "Error, no *_census_stage_0.ini file found in %s", self.meteor.fastq_dir
            )
            sys.exit()
        else:
            # Not sure if it's a good idea to delete a temporary file
            rmtree(self.meteor.tmp_dir, ignore_errors=True)
        return True
