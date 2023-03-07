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
from dataclasses import dataclass, field
from tempfile import mkdtemp, NamedTemporaryFile
import tempfile
from configparser import ConfigParser
from pathlib import Path
from subprocess import check_call
from meteor.mapper import Mapper
from meteor.mapper2 import Mapper2
from meteor.session import Session, Component
from typing import Type
from collections import defaultdict
import itertools
from pysam import index, idxstats, AlignmentFile
from time import perf_counter


@dataclass
class Counter(Session):
    """Counter session map and count
    """
    meteor: Type[Component]
    counting_type: str
    mapping_type: str
    trim: int
    alignment_number: int
    counting_only: bool
    mapping_only: bool
    pysam_test: bool = True
    ini_data: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.meteor.tmp_dir = Path(mkdtemp(dir=self.meteor.tmp_path))
        self.meteor.mapping_dir.mkdir(exist_ok=True)

    def count_index_fastq(self, fastq_file: Path, output_desc: tempfile._TemporaryFileWrapper) -> tuple:
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
            for read_count, line in enumerate(in_fq, start=1):  # pylint: disable=unused-variable
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
        """Write configuration file for reference genome

        :return: (ConfigParser) A configparser object
        """
        config = ConfigParser()
        config["worksession"] = {
            "meteor.reference.dir": self.meteor.ref_dir.name,
            "meteor.db.type": "binary",
            "meteor.mapping.program": "bowtie2",
            "meteor.mapping.file.format": "sam",
            "meteor.is.cpu.percentage": "0",
            "meteor.cpu.count": str(self.meteor.threads),
            "meteor.excluded.reference.count": "0"
        }
        config["main_reference"] = {
            "meteor.reference.name": self.meteor.ref_name,
            "meteor.matches": str(self.alignment_number),
            "meteor.mismatches": "5",
            "meteor.is.perc.mismatches": "1",
            "meteor.bestalignment": "1",
            "meteor.mapping.prefix.name": f"mapping_vs_{ref_ini['reference_info']['reference_name']}",
            "meteor.counting.prefix.name": f"vs_{ref_ini['reference_info']['reference_name']}"
        }
        return config

    def launch_mapping(self) -> None:
        """Create temporary indexed files and map aga
        """
        logging.info("Launch mapping")
        # loop on each library
        for library, dict_data in self.ini_data.items():
            census = dict_data["census"]
            sample_info = census["sample_info"]  # reference
            sample_file = census["sample_file"]
            logging.info("Meteor Mapping task description")
            logging.info("Sample name = %s",  sample_info["sample_name"])
            logging.info("Library name = %s", sample_info["full_sample_name"])
            logging.info("Project name = %s", sample_info["project_name"])
            logging.info("Sequencing device = %s", sample_info["sequencing_device"])
            logging.info("Workflow = %s", library.name)

            # reindexing this library reads and fill FLibraryIndexerReport
            fastq_path = self.meteor.fastq_dir / sample_file["fastq_file"]
            try:
                with NamedTemporaryFile(mode="wt", dir=self.meteor.tmp_dir) as output_desc:
                    read_count, base_count = self.count_index_fastq(fastq_path, output_desc)
                    census.set("sample_info", "census_status",  str(1))
                    census.set("sample_info", "indexed_read_length", str(1))
                    census.set("sample_info", "sequenced_read_count",  str(read_count))
                    census.set("sample_info", "indexed_sequenced_read_count",  str(read_count))
                    census.set("sample_info", "indexed_sequenced_base_count", str(base_count))
                    census.set("sample_info", "is_data_prefixed", str(0))
                    dict_data["census"] = census
                    # mapping this library on the reference
                    mapping_process = Mapper(self.meteor, dict_data,
                                             Path(output_desc.name),
                                             self.mapping_type, self.trim,
                                             self.alignment_number)
                    if not mapping_process.execute():
                        raise ValueError(f"Error, TaskMainMapping failed: {library}")
            except IOError:
                logging.error("Cannot create temporary files in %s", self.meteor.tmp_dir.name)

    def launch_counting(self, workflow_ini: Path) -> None:
        """Launch meteor counter

        :param workflow_ini: A path object to workflow configuration
        """
        logging.info("Launch counting")
        # "-t", str(self.meteor.tmp_dir) + "/"
        start = perf_counter()
        check_call(["meteor-counter", "-i", str(self.meteor.fastq_dir) + "/",
                    "-p", str(self.meteor.ref_dir.parent.resolve()) + "/",
                    "-o", str(self.meteor.mapping_dir) + "/", "-w", workflow_ini,
                    "-c", self.counting_type, "-f"])
        print(f"Meteor counter completed Execution in {perf_counter() - start} seconds")

    def launch_mapping2(self) -> None:
        """Create temporary indexed files and map aga
        """
        logging.info("Launch mapping")
        list_fastq_path = []
        # loop on each library
        for library, dict_data in self.ini_data.items():  # pylint: disable=unused-variable
            census = dict_data["census"]
            sample_file = census["sample_file"]
            # reindexing this library reads and fill FLibraryIndexerReport
            list_fastq_path += [str(self.meteor.fastq_dir / sample_file["fastq_file"])]
        else:
            # mapping this library on the reference
            mapping_process = Mapper2(self.meteor, dict_data, list_fastq_path,
                                      self.mapping_type, self.trim,
                                      self.alignment_number, self.counting_type)
            if not mapping_process.execute():
                raise ValueError("Error, TaskMainMapping failed")

    def write_table(self, bamfile: str, outfile: Path) -> None:
        """Function that create a count table using pysam. First index the BAM file,
        then count reads using the function idxstats from pysam, and output a count
        table.

        :param bamfile: (Path) BAM file to count
        :outfile: (str) count table name
        """
        # index the bam file
        index(bamfile)
        # create count table
        table = idxstats(bamfile)
        # write the count table
        with outfile.open("wt", encoding="UTF-8") as out:
            for line in table:
                out.write(line)

    def filter_bam(self, bamdesc):
        """Function that count reads from a BAM file, using the given methods in count:
        "total" or "shared". If bam is set to 'True', a BAM file containing only
        the used alignment for the counting is generated.

        :param bamfile [STR] = BAM file to count

        Returns:
            if count = total : new_filename [STR] = BAM file
            if count = shared :
                    database [DICT] = contains length of reference genomes.
                                        key : reference genome
                                        value : size
                    tmp_genomes [DICT] =
        """
        tmp_score = {}
        genomes = defaultdict(list)
        # contains a list a alignment of each read
        reads = defaultdict(list)
        for element in bamdesc:
            # if not element.is_unmapped: # if read mapped
            if not element.has_tag("AS"):
                raise ValueError("Missing 'AS' field.")
            if element.is_read1:
                # get read orientation
                read_id = f"{element.qname}_1"
            else:
                read_id = f"{element.qname}_2"
            # get alignment score
            score = element.get_tag("AS")
            # get previous score for the read
            prev_score = tmp_score.get(read_id, score)

            # if same score
            if prev_score == score:
                tmp_score[read_id] = score
                # add the genome to the list if it doesn't exist
                if element.reference_name not in genomes[read_id]:
                    genomes[read_id].append(element.reference_name)
                    if self.counting_type == "total_reads":
                        # append line for filtered BAM file
                        reads[read_id].append(element)
            # if previous score is lower
            elif prev_score < score:
                # set the new score
                tmp_score[read_id] = score
                # reinitialise all
                genomes[read_id] = [element.reference_name]
                if self.counting_type == "total_reads":
                    reads[read_id] = [element]
        return reads, genomes

    def uniq_from_mult(self, genome_dict, database):
        """
        Function that filter unique reads from all reads. Multiple reads are
        reads that map to more than one genome. And Unique reads are reads that map
        only on one genome.

        Args:
            genome_dict [DICT] = dictionary containing reads as key and a list of
                                genome where read mapped as value
            unique_dict [DICT] = contains for each reference genome the number of
                                unique reads

        Returns:
            genome_dict [DICT] = the dictionary without unique reads
            unique_dict [DICT] = nb of unique read of each reference genome
        """
        unique_reads = []
        unique_dict = dict.fromkeys(database, 0)
        for key in genome_dict:
            # if read mapp with best score to only one genome:
            if len(genome_dict[key]) != 1:
                continue
            # get genome id
            genome = "".join(genome_dict[key])
            # add to unique read dictionnary
            unique_dict[genome] += 1
            unique_reads.append(key)

        for key in unique_reads:
            # delete the key from TMP dictionnary
            del genome_dict[key]

        return genome_dict, unique_dict

    def compute_co(self, genome_dict, unique_dict):
        """
        Compute genome specific coefficient "Co" for each multiple read.

        Args:
            genome_dict [DICT] = Contains for each multiple read, all the reference
                                genomes name where he mapped
            unique_dict [DICT] = nb of unique read of each reference genome

        Returns:
            Co [DICT] = contains coefficient values for each couple (read, genome)
            read_dict [DICT] = contains all the multiple reads of a genome
                                key : genome
                                value : list a multiple reads
        """
        co = {}
        read_dict = {}

        for key in genome_dict:
            # for a multiple read, gets nb of unique reads from each genome
            s = [unique_dict[genome] for genome in genome_dict[key]]
            # total number of unique reads
            som = reduce(lambda x, y: x+y, s)
            if som == 0:
                continue
            for genome in genome_dict[key]:
                # get the nb of unique reads
                nb_unique = unique_dict[genome]
                if nb_unique == 0:
                    continue
                # calculate Co of multiple read for the given genome
                co[(key, genome)] = nb_unique / float(som)
                # append to the dict
                read_dict.setdefault(genome, []).append(key)

        # get all the multiple reads of a genome
        for genome, reads in read_dict.items():
            read_dict[genome] = list(set(reads))

        return read_dict, co

    def get_co_coefficient(self, genome, read_dict, co_dict):
        for read in read_dict[genome]:
            yield co_dict[(read, genome)]

    def compute_abm(self, read_dict, co_dict, database):
        """
        Compute multiple reads abundance for each genome.
        Multiple reads abundance of a genome is equal to the sum of all Co
        coefficient

        Args:
            read_dict [DICT] = list of multiple read for each reference genome
            co_dict [DICT] = contains coefficient values for each couple
                            (read, genome)
            multiple_dict [DICT] = abundance of multiple reads for each genome

        Returns:
            multiple_dict [DICT] = abundance of multiple reads for each genome
        """
        multiple_dict = dict.fromkeys(database, 0)
        for genome in read_dict:
            # get a list of all the Co coefficient for each reference genome
            # sum them, (we obtain the abundance)
            multiple_dict[genome] = sum(self.get_co_coefficient(genome, read_dict, co_dict))
        return multiple_dict

    def compute_abs(self, unique_dict, multiple_dict):
        """
        Compute the abundance of a each genome.
            Abundance = Abundance unique reads + Abundance multiple reads

        Args:
            unique_dict [DICT] = nb of unique read of each reference genome
            multiple_dict [DICT] = abundance of multiple reads for each genome

        Returns:
            abundance_dict [DICT] = contains abundance of each reference genome
        """
        return {genome: unique_dict[genome] + multiple_dict[genome] for genome in unique_dict}

    def write_stat(self, output: Path, abundance_dict: dict, database):
        """
        Write count table.

        Args:
            output [STRING] = output filename
            abundance_dict [DICT] = contains abundance of each reference genome
            database [DICT] = contrains length of each reference genome

        No Returns
        """
        with open(output, "wt", encoding="UTF-8") as out:
            for genome, abundance in abundance_dict.items():
                out.write(f"{genome}\t{database[genome]}\t{abundance}\n")

    def launch_counting2(self, bam_file: Path, count_file: Path) -> None:
        """Launch meteor counter

        :param bam_file: (Path) A path to the input bam file
        :param count_file: (Path) A path to the output count file
        """
        logging.info("Launch counting")
        # "-t", str(self.meteor.tmp_dir) + "/"
        # open the BAM file
        with AlignmentFile(str(bam_file.resolve()), "rb") as bamdesc:
            if self.mapping_type == "total_reads":
                # name of the filtered BAM file
                new_filename = bam_file.parent / f"unsorted_filtered_{bam_file.name}"
                reads, genomes = self.filter_bam(bamdesc)
                read_list = list(itertools.chain(reads.values()))
                merged_list = list(itertools.chain.from_iterable(read_list))
                # writing the filtered BAM file
                with AlignmentFile(str(new_filename.resolve()), "wb", template=bamdesc) as total_reads:
                    for element in merged_list:
                        total_reads.write(element)
                self.write_table(str(new_filename.resolve()), count_file)
            elif self.counting_type == "smart_shared_reads":
                # create a dictionary containing the length of reference genomes
                # get name of reference sequence
                references = bamdesc.references
                # get reference length
                lengths = bamdesc.lengths
                database = dict(zip(references, lengths))
                reads, genomes = self.filter_bam(bamdesc)
                # writing a bam file for total_reads counting and if a bam=True
                # returns reads and alignements informations if shared counting
                for key, values in genomes.items():
                    genomes[key] = list(set(values))
                # Filter unique reads from multiple reads
                genomes, unique = self.uniq_from_mult(genomes, database)
                # For multiple reads compute Co
                reads, coef_read = self.compute_co(genomes, unique)
                # calculate abundance
                # TODO unique reads ignored here
                multiple = self.compute_abm(reads, coef_read, database)
                # Calculate reference abundance & write count table
                abundance = self.compute_abs(unique, multiple)
                self.write_stat(count_file, abundance, database)

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
            logging.error("Error, no *_reference.ini file found in %s. "
                          "One *_reference.ini is expected", self.meteor.ref_dir)
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
                    if self.pysam_test:
                        stage1_dir = self.meteor.mapping_dir / sample_info["sample_name"]
                    else:
                        stage1_dir = (
                            self.meteor.mapping_dir / sample_info["sample_name"] /
                            f"mapping_vs_{self.meteor.ref_name}_{census_ini['sample_info']['full_sample_name']}"
                        )
                    # sample_info["full_sample_name"]
                    stage1_dir.mkdir(exist_ok=True, parents=True)
                    self.ini_data[library] = {
                        "census": census_ini,
                        "directory": stage1_dir,
                        "Stage1FileName": stage1_dir / library.name.replace("stage_0", "stage_1"),
                        "reference": ref_ini
                    }
                if not self.ini_data[library]["Stage1FileName"].exists():
                    mapping_done = False
            else:
                if not self.counting_only:
                    # mapping already done and no overwriting
                    if mapping_done:
                        logging.info("Mapping already done for sample: %s", sample_info["sample_name"])
                        logging.info("Skipped !")
                    else:
                        if self.pysam_test:
                            self.launch_mapping2()
                        else:
                            self.launch_mapping()
                if not self.mapping_only:
                    logging.info("Launch mapping")
                    config = self.set_workflow_config(ref_ini)
                    workflow_ini = self.meteor.mapping_dir / "workflow.ini"
                    self.save_config(config, workflow_ini)
                    bam_file = self.ini_data[library]["directory"] / f"{sample_info['sample_name']}.bam"
                    count_file = self.ini_data[library]["directory"] / f"{sample_info['sample_name']}.tsv"
                    if self.pysam_test:
                        start = perf_counter()
                        if self.counting_type == "unique_reads":
                            self.write_table(str(bam_file.resolve()), count_file)
                        else:
                            self.launch_counting2(bam_file, count_file)
                        logging.info("Completed counting creation in %f seconds", perf_counter() - start)
                    else:
                        self.launch_counting(workflow_ini)
            logging.info("Done ! Job finished without errors ...")
            self.meteor.tmp_dir.rmdir()
        except AssertionError:
            logging.error("Error, no *_census_stage_0.ini file found in %s", self.meteor.fastq_dir)
            sys.exit()
        return True
