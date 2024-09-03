#!/bin/env python3
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

"""Meteor - A plateform for quantitative metagenomic profiling of complex ecosystems"""


import sys
import logging
from argparse import ArgumentParser, ArgumentTypeError, Namespace, RawTextHelpFormatter
from pathlib import Path
from meteor.session import Component
from meteor.fastqimporter import FastqImporter
from meteor.referencebuilder import ReferenceBuilder
from meteor.mapper import Mapper
from meteor.counter import Counter
from meteor.downloader import Downloader
from meteor.profiler import Profiler
from meteor.merging import Merging
from meteor.strain import Strain
from meteor.treebuilder import TreeBuilder
from tempfile import TemporaryDirectory
from importlib.metadata import version


class Color:
    BOLD = "\033[1m"
    END = "\033[0m"


def get_logging() -> logging.Logger:  # pragma: no cover
    """Configure a stream and file logging

    :return: (logging.logger) A logger object
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s :: %(levelname)s :: %(message)s")
    # Create log file
    file_handler = logging.FileHandler("meteor.log", "a", encoding="UTF-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Stream in the the console
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    return logger


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if not myfile.exists():
            msg = f"{myfile.name} does not exist."
        else:
            msg = f"{myfile.name} is not a regular file."
        raise ArgumentTypeError(msg)
    return myfile


def isdir(path: str) -> Path:  # pragma: no cover
    """Check if path can be valid directory.

    :param path: Path to the directory

    :raises ArgumentTypeError: If directory does not exist

    :return: (str) Path object of the directory
    """
    mydir = Path(path)
    if mydir.is_file():
        msg = f"{mydir.name} is a file."
        raise ArgumentTypeError(msg)
    return mydir


def isborned01(x: str) -> float:
    """Check if a float is comprised between 0 and 1.

    :param x: float number to check
    :raises ArgumentTypeError: If x is > 1 or < 0
    :return: (float) float number
    """
    # print(x)
    # if not isinstance(x, float):
    #     raise ArgumentTypeError("Value must be a numerical")
    x_float = float(x)
    if x_float < 0.0 or x_float > 1.0:
        msg = "Should be comprised between 0 and 1."
        raise ArgumentTypeError(msg)
    return x_float


def num_threads(value):

    try:
        value = int(value)
    except ValueError as value_err:
        raise ArgumentTypeError(
            "the number of threads is not an integer"
        ) from value_err

    if value <= 0:
        raise ArgumentTypeError("the minimum number of threads is 1")
    return value


def get_arguments() -> Namespace:  # pragma: no cover
    """
    Meteor help and arguments

    No arguments

    Return : parser
    """
    parser = ArgumentParser(description=Color.BOLD + __doc__ + Color.END, prog="Meteor")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s version { version('meteor')}"
    )
    subparsers = parser.add_subparsers(
        title="positional arguments",
        help="Select activity",
        dest="command",
        required=True,
    )
    # Define download argument parsing
    download_parser = subparsers.add_parser("download", help="Download a catalogue")
    download_parser.add_argument(
        "-i",
        dest="user_choice",
        type=str,
        required=True,
        choices=Downloader.get_available_catalogues(),
        help="Select the catalogue to download.",
    )
    download_parser.add_argument(
        "--fast",
        dest="taxonomy",
        action="store_true",
        help="Select the short catalogue variant (only for taxonomic profiling).",
    )
    download_parser.add_argument(
        "-o",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Directory where the downloaded catalogue is saved.",
    )
    download_parser.add_argument(
        "-c",
        dest="check_md5",
        action="store_true",
        help="Check the md5sum of the catalogue after download.",
    )
    # Define download argument parsing
    reference_parser = subparsers.add_parser(
        "build", help="Import a custom gene catalogue"
    )
    reference_parser.add_argument(
        "-i",
        dest="input_fasta_file",
        type=isfile,
        required=True,
        help="Fasta file of the custom gene catalogue (gzip, bzip2 and xz compression accepted).",
    )
    reference_parser.add_argument(
        "-o",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Directory where the catalogue and its index are saved.",
    )
    reference_parser.add_argument(
        "-n",
        dest="ref_name",
        metavar="REFERENCE_NAME",
        type=str,
        required=True,
        help="Name of the gene catalogue (ansi-string without space).",
    )
    reference_parser.add_argument(
        "-t",
        dest="threads",
        default=ReferenceBuilder.DEFAULT_NUM_THREADS,
        type=num_threads,
        help="Number of threads to launch while indexing the catalogue (default: %(default)d).",
    )
    # Define fastq argument parsing
    fastq_parser = subparsers.add_parser(
        "fastq", formatter_class=RawTextHelpFormatter, help="Import fastq files"
    )
    fastq_parser.add_argument(
        "-i",
        dest="input_fastq_dir",
        type=isdir,
        required=True,
        help="Directory containing all input fastq files with .fastq or .fq. extensions "
        "(gzip, bzip2 and xz compression accepted).\n"
        "Paired-end files must be named : file_R1.[fastq/fq] & file_R2.[fastq/fq] "
        "or file_1.[fastq/fq] & file_2.[fastq/fq]",
    )
    fastq_parser.add_argument(
        "-p",
        dest="ispaired",
        action="store_true",
        help="Fastq files are paired.",
    )
    fastq_parser.add_argument(
        "-m",
        dest="mask_sample_name",
        type=str,
        help="Regular expression (between quotes) for extracting sample name.",
    )
    fastq_parser.add_argument(
        "-o",
        dest="fastq_dir",
        type=isdir,
        required=True,
        help="Directory where the fastq repository is created.",
    )
    # Define mapping argument parsing
    mapping_parser = subparsers.add_parser(
        "mapping",
        formatter_class=RawTextHelpFormatter,
        help="Map reads against a gene catalogue and compute raw gene counts",
    )
    mapping_parser.add_argument(
        "-i",
        dest="fastq_dir",
        type=isdir,
        required=True,
        help="Directory corresponding to the sample to process.\n"
        "(contains sequencing metadata files ending with _census_stage_0.json)",
    )
    mapping_parser.add_argument(
        "-r",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Directory corresponding to the gene catalog against which reads are mapped.\n"
        "(contains a file ending with *_reference.json)",
    )
    mapping_parser.add_argument(
        "-o",
        dest="mapping_dir",
        type=isdir,
        required=True,
        help="Directory where mapping and raw gene counts of the sample are saved.",
    )
    mapping_parser.add_argument(
        "-p",
        dest="mapping_type",
        type=str,
        choices=Mapper.MAPPING_TYPES,
        default=Mapper.DEFAULT_MAPPING_TYPE,
        help="Strategy to map reads against the catalogue (default: %(default)s).",
    )
    mapping_parser.add_argument(
        "--trim",
        dest="trim",
        type=int,
        default=Mapper.DEFAULT_TRIM,
        help="Trim reads exceeding TRIM bases before mapping "
        f"(default: %(default)d).\nIf {Mapper.NO_TRIM}, no trim.",
    )
    mapping_parser.add_argument(
        "--align",
        dest="alignment_number",
        type=int,
        default=Mapper.DEFAULT_ALIGNMENT_NUMBER,
        help="Maximum number alignments to report for each read (default: %(default)d)",
    )
    mapping_parser.add_argument(
        "-c",
        dest="counting_type",
        type=str,
        default=Counter.DEFAULT_COUNTING_TYPE,
        choices=Counter.COUNTING_TYPES,
        help="Strategy to calculate raw gene counts (default: %(default)s).",
    )
    mapping_parser.add_argument(
        "--core_size",
        dest="core_size",
        type=int,
        default=Component.DEFAULT_CORE_SIZE,
        help="Number of signature genes per species (MSP) used to estimate their respective abundance "
        "(default: %(default)d).",
    )
    mapping_parser.add_argument(
        "--id",
        dest="identity_threshold",
        type=isborned01,
        # default=Counter.DEFAULT_IDENTITY_THRESHOLD,
        help="Select only read alignments with a nucleotide identity >= IDENTITY_THRESHOLD "
        f"(default: auto).\nIf {Counter.NO_IDENTITY_THRESHOLD}, no filtering.",
    )
    mapping_parser.add_argument(
        "--ka",
        dest="keep_all_alignments",
        action="store_true",
        help="Keep raw bowtie2 output in cram format. "
        "Required for calculating gene counts with another strategy.",
    )
    mapping_parser.add_argument(
        "--kf",
        dest="keep_filtered_alignments",
        action="store_true",
        help="Keep filtered alignments on marker genes. Required for strain analysis.",
    )
    mapping_parser.add_argument(
        "--tmp",
        dest="tmp_path",
        type=isdir,
        help="Directory where temporary files (e.g. cram) are stored",
    )
    mapping_parser.add_argument(
        "-t",
        dest="threads",
        type=num_threads,
        default=Mapper.DEFAULT_NUM_THREADS,
        help="Number of alignment threads to launch (default: %(default)d).",
    )
    # Define profiler argument parsing
    profiling_parser = subparsers.add_parser(
        "profile",
        formatter_class=RawTextHelpFormatter,
        help="Compute species and functional abundance tables",
    )
    profiling_parser.add_argument(
        "-i",
        dest="mapped_sample_dir",
        required=True,
        type=isdir,
        help="Directory with raw gene counts of the sample to process.\n"
        "(contains a metadata file ending with _census_stage_1.json)",
    )
    profiling_parser.add_argument(
        "-r",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Directory corresponding to the catalog used to generate raw gene counts.\n"
        "(contains a file ending with *_reference.json)",
    )
    profiling_parser.add_argument(
        "-o",
        dest="profiled_sample_dir",
        type=isdir,
        required=True,
        help="Directory where species and functional abundance tables of the sample are saved.",
    )
    profiling_parser.add_argument(
        "-l",
        dest="rarefaction_level",
        type=int,
        default=Profiler.DEFAULT_RAREFACTION_LEVEL,
        help=f"Rarefaction level. If <= {Profiler.NO_RAREFACTION}, no rarefation is performed "
        "(default: %(default)d).",
    )
    profiling_parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=Profiler.DEFAULT_RANDOM_SEED,
        help="Seed of the random number generator used for rarefaction "
        "(default: %(default)d).",
    )
    profiling_parser.add_argument(
        "-n",
        dest="normalization",
        type=str,
        choices=Profiler.NORMALIZATIONS,
        default=Profiler.DEFAULT_NORMALIZATION,
        help="Normalization applied to raw gene counts (default: %(default)s).",
    )
    profiling_parser.add_argument(
        "-c",
        dest="coverage_factor",
        type=float,
        default=Profiler.DEFAULT_COVERAGE_FACTOR,
        help="Multiplication factor for coverage normalization "
        "(default: %(default).0f).",
    )
    profiling_parser.add_argument(
        "--core_size",
        dest="core_size",
        type=int,
        default=Component.DEFAULT_CORE_SIZE,
        help="Number of signature genes per species (MSP) used to estimate their respective abundance "
        "(default: %(default)d).",
    )
    profiling_parser.add_argument(
        "--msp_filter",
        dest="msp_filter",
        type=isborned01,
        # default=Profiler.DEFAULT_MSP_FILTER,
        help="Minimal proportion of core genes detected in a sample to consider a species (MSP) as present "
        "(default: auto).",
    )
    profiling_parser.add_argument(
        "--completeness",
        type=isborned01,
        default=Profiler.DEFAULT_COMPLETENESS,
        help="Cutoff above which a module is considered as present in a species.\n"
        "Value between 0.0 and 1.0 (default: %(default).1f)."
        "",
    )
    # Define merging argument parsing
    merging_parser = subparsers.add_parser(
        "merge",
        formatter_class=RawTextHelpFormatter,
        help="Merge gene, species and functional abundance tables from multiple samples",
    )
    merging_parser.add_argument(
        "-i",
        dest="profile_dir",
        required=True,
        type=isdir,
        help="Directory containing subdirectories (one per sample) with abundance tables to be merged.\n"
        "(each subdirectory contains a metadata file ending with _census_stage_2.json)",
    )
    merging_parser.add_argument(
        "-r",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Directory corresponding to the gene catalog used to generate the abundance tables.\n"
        "(contains a file ending with *_reference.json)",
    )
    merging_parser.add_argument(
        "-a",
        dest="min_msp_abundance",
        type=float,
        default=Merging.DEFAULT_MIN_MSP_ABUNDANCE,
        help="Minimum msp abundance (default >= %(default).1f)",
    )
    merging_parser.add_argument(
        "-n",
        dest="min_msp_occurrence",
        type=int,
        default=Merging.DEFAULT_MIN_MSP_OCCURRENCE,
        help="Report only species (MSPs) occuring in at least n samples "
        "(default: %(default)d).",
    )
    merging_parser.add_argument(
        "-s",
        dest="remove_sample_with_no_msp",
        action="store_true",
        help="Remove samples with no detected species (MSPs) "
        "(default: %(default)s).",
    )
    # merging_parser.add_argument(
    #    "-m",
    #    dest="output_mpa",
    #    action="store_true",
    #    help="Save the merged species abundance table in the style of MetaPhlan "
    #    "(default: %(default)s).",
    # )
    merging_parser.add_argument(
        "-b",
        dest="output_biom",
        action="store_true",
        help="Save the merged species abundance table in biom format "
        "(default: %(default)s).",
    )
    # merging_parser.add_argument(
    #    "--tax_lev",
    #    dest="taxonomic_level",
    #    default=Merging.DEFAULT_MPA_TAXONOMIC_LEVEL,
    #    choices=Merging.MPA_TAXONOMIC_LEVELS,
    #    help="""The taxonomic level for mpa output (default: %(default)s):
    #                    'a' : all taxonomic levels
    #                    'k' : kingdoms
    #                    'p' : phyla only
    #                    'c' : classes only
    #                    'o' : orders only
    #                    'f' : families only
    #                    'g' : genera only
    #                    's' : species only
    #                    't' : MSPs only""",
    # )
    merging_parser.add_argument(
        "-o",
        dest="merging_dir",
        required=True,
        type=isdir,
        help="Directory where the merged abundance tables are saved.",
    )
    merging_parser.add_argument(
        "-p",
        dest="prefix",
        default=Merging.DEFAULT_PREFIX,
        help="Prefix added to output filenames " "(default: %(default)s).",
    )
    merging_parser.add_argument(
        "-g",
        dest="output_gene_matrix",
        action="store_true",
        help="Merge gene abundance tables.",
    )
    strain_parser = subparsers.add_parser(
        "strain", help="Identifies strains from metagenomic samples"
    )
    strain_parser.add_argument(
        "-i",
        dest="mapped_sample_dir",
        required=True,
        type=isdir,
        help="Path to the mapped sample directory.",
    )
    strain_parser.add_argument(
        "-r",
        dest="ref_dir",
        type=isdir,
        required=True,
        help="Path to reference directory (Path containing *_reference.json)",
    )
    strain_parser.add_argument(
        "-d",
        dest="max_depth",
        default=Strain.DEFAULT_MAX_DEPTH,
        type=int,
        help="Maximum number of reads taken in account (default: %(default)d).",
    )
    strain_parser.add_argument(
        "-p",
        dest="min_depth",
        default=Strain.DEFAULT_MIN_DEPTH,
        type=int,
        help=f"""Minimum depth (default: >= %(default)d).
        Values should be comprised between {Strain.MIN_DEPTH} and the maximum depth
        ({Strain.MAX_MIN_SNP_DEPTH} reads are taken in account).""",
    )
    strain_parser.add_argument(
        "-s",
        dest="min_snp_depth",
        default=Strain.DEFAULT_MIN_SNP_DEPTH,
        choices=range(Strain.MIN_MIN_SNP_DEPTH, Strain.MAX_MIN_SNP_DEPTH + 1),
        metavar="MIN_SNP_DEPTH",
        type=int,
        help=f"""Minimum snp depth (default: >= %(default)d).
        Values should be comprised between {Strain.MIN_MIN_SNP_DEPTH} and the maximum depth
        ({Strain.MAX_MIN_SNP_DEPTH} reads are taken in account).""",
    )
    strain_parser.add_argument(
        "-f",
        dest="min_frequency",
        default=Strain.DEFAULT_MIN_FREQUENCY,
        type=isborned01,
        help="Minimum frequency for alleles (default: >= %(default).2f).",
    )
    strain_parser.add_argument(
        "-l",
        dest="ploidy",
        default=Strain.DEFAULT_PLOIDY,
        type=int,
        help="Ploidy (default: %(default)d).",
    )
    strain_parser.add_argument(
        "-m",
        dest="min_msp_coverage",
        default=Strain.DEFAULT_MIN_MSP_COVERAGE,
        choices=range(Strain.MIN_MIN_MSP_COVERAGE, Strain.MAX_MIN_MSP_COVERAGE + 1),
        metavar="MIN_MSP_COVERAGE",
        type=int,
        help=f"""Minimum number of genes from the MSP that are covered (default: >= %(default)d).
        Values should be comprised between {Strain.MIN_MIN_MSP_COVERAGE} and {Strain.MAX_MIN_MSP_COVERAGE}""",
    )
    strain_parser.add_argument(
        "-c",
        dest="min_gene_coverage",
        default=Strain.DEFAULT_MIN_GENE_COVERAGE,
        type=isborned01,
        help="Minimum gene coverage from 0 to 1 (default: >= %(default).1f).",
    )
    strain_parser.add_argument(
        "--core_size",
        dest="core_size",
        type=int,
        default=Component.DEFAULT_CORE_SIZE,
        help="Number of signature genes per species (MSP) used to estimate their respective abundance "
        "(default: %(default)d).",
    )
    strain_parser.add_argument(
        "-o",
        dest="strain_dir",
        type=isdir,
        required=True,
        help="Path to output directory.",
    )
    strain_parser.add_argument(
        "--kc",
        dest="keep_consensus",
        action="store_true",
        help="Keep consensus marker genes (default: False, set to True to recompute strain)",
    )
    strain_parser.add_argument(
        "--tmp",
        dest="tmp_path",
        type=isdir,
        help="Path to the directory where temporary files are stored",
    )
    strain_parser.add_argument(
        "-t",
        dest="threads",
        default=Strain.DEFAULT_NUM_THREADS,
        type=num_threads,
        help="Number of threads to perform variant calling (default: %(default)d).",
    )
    tree_parser = subparsers.add_parser(
        "tree", help="Compute phylogenetical tree from detected strains"
    )
    tree_parser.add_argument(
        "-i",
        dest="strain_dir",
        required=True,
        type=isdir,
        help="Path to the strain directory.",
    )
    tree_parser.add_argument(
        "-g",
        dest="max_gap",
        default=TreeBuilder.DEFAULT_MAX_GAP,
        type=isborned01,
        help="Removes sites constitued of >= cutoff gap character (default: >= %(default).1f).",
    )
    tree_parser.add_argument(
        "-s",
        dest="min_info_sites",
        default=TreeBuilder.DEFAULT_MIN_INFO_SITES,
        type=isborned01,
        help="Minimum number of informative sites in the alignment (default: >= %(default)d).",
    )
    tree_parser.add_argument(
        "-f",
        dest="format",
        default=TreeBuilder.DEFAULT_OUTPUT_FORMAT,
        choices=TreeBuilder.OUTPUT_FORMATS,
        type=str,
        help="Output image format (default: %(default)s).",
    )
    tree_parser.add_argument(
        "-w",
        dest="width",
        default=TreeBuilder.DEFAULT_WIDTH,
        type=int,
        help="Output image width (default: %(default)dpx).",
    )
    tree_parser.add_argument(
        "-H",
        dest="height",
        default=TreeBuilder.DEFAULT_HEIGHT,
        type=int,
        help="Output image height (default: %(default)dpx).",
    )
    tree_parser.add_argument(
        "-o",
        dest="tree_dir",
        type=isdir,
        required=True,
        help="Path to output directory.",
    )
    tree_parser.add_argument(
        "--tmp",
        dest="tmp_path",
        type=isdir,
        help="Path to the directory where temporary files are stored",
    )
    tree_parser.add_argument(
        "-t",
        dest="threads",
        default=TreeBuilder.DEFAULT_NUM_THREADS,
        type=num_threads,
        help="Number of threads when infering each tree (default: %(default)d).",
    )
    subparsers.add_parser("test", help="Test meteor installation")
    return parser.parse_args(args=None if sys.argv[1:] else ["--help"])


def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Let us logging
    logger = get_logging()
    # Get arguments
    args = get_arguments()
    # version = importlib.metadata.version("meteor")
    # print("Meteor version", version)
    # Create a meteor dataset
    meteor = Component
    # Import FASTQ
    if args.command == "fastq":
        meteor.fastq_dir = args.fastq_dir
        # args.isdispatched,
        fastq_importer = FastqImporter(
            meteor,
            args.input_fastq_dir,
            args.ispaired,
            args.mask_sample_name,
        )
        fastq_importer.execute()
    # Import reference
    elif args.command == "build":
        meteor.ref_name = args.ref_name
        meteor.ref_dir = args.ref_dir
        meteor.threads = args.threads
        # pysam_test=args.pysam_test
        reference_builder = ReferenceBuilder(
            meteor=meteor, input_fasta=args.input_fasta_file
        )
        reference_builder.execute()
    # Run mapping
    elif args.command == "mapping":
        meteor.fastq_dir = args.fastq_dir
        meteor.mapping_dir = args.mapping_dir
        meteor.ref_dir = args.ref_dir
        meteor.tmp_path = args.tmp_path
        meteor.threads = args.threads
        # args.pysam_test
        counter = Counter(
            meteor,
            args.counting_type,
            args.mapping_type,
            args.trim,
            args.identity_threshold,
            args.alignment_number,
            args.core_size,
            args.keep_all_alignments,
            args.keep_filtered_alignments,
        )
        counter.execute()
    # Run strain
    elif args.command == "strain":
        meteor.mapped_sample_dir = args.mapped_sample_dir
        meteor.ref_dir = args.ref_dir
        meteor.tmp_path = args.tmp_path
        meteor.threads = args.threads
        meteor.strain_dir = args.strain_dir
        strain_id = Strain(
            meteor,
            args.max_depth,
            args.min_depth,
            args.min_snp_depth,
            args.min_frequency,
            args.ploidy,
            args.min_msp_coverage,
            args.min_gene_coverage,
            args.core_size,
            args.keep_consensus,
        )
        strain_id.execute()
    # Compute trees
    elif args.command == "tree":
        meteor.strain_dir = args.strain_dir
        meteor.tree_dir = args.tree_dir
        meteor.threads = args.threads
        meteor.tmp_path = args.tmp_path
        trees = TreeBuilder(
            meteor,
            args.max_gap,
            args.min_info_sites,
            args.width,
            args.height,
            args.format,
        )
        trees.execute()
    # Run download catalogues
    elif args.command == "download":
        meteor.ref_name = args.user_choice
        meteor.ref_dir = args.ref_dir
        downloader = Downloader(meteor, args.user_choice, args.taxonomy, args.check_md5)
        downloader.execute()
    # Run profiling
    elif args.command == "profile":
        meteor.mapping_dir = args.mapped_sample_dir
        meteor.profile_dir = args.profiled_sample_dir
        meteor.ref_dir = args.ref_dir
        profiler = Profiler(
            meteor,
            args.rarefaction_level,
            args.seed,
            args.normalization,
            args.core_size,
            args.msp_filter,
            args.completeness,
            args.coverage_factor,
        )
        profiler.execute()
    # Run merging
    elif args.command == "merge":
        meteor.ref_dir = args.ref_dir
        meteor.profile_dir = args.profile_dir
        meteor.merging_dir = args.merging_dir
        merging = Merging(
            meteor,
            args.prefix,
            args.min_msp_abundance,
            args.min_msp_occurrence,
            args.remove_sample_with_no_msp,
            False,
            None,
            # args.output_mpa,
            # args.taxonomic_level,
            args.output_biom,
            args.output_gene_matrix,
        )
        merging.execute()
    # Testing
    else:
        with TemporaryDirectory() as tmpdirname:
            meteor.ref_name = "test"
            meteor.ref_dir = Path(tmpdirname) / "ref"
            meteor.threads = 1
            meteor.tmp_path = Path(tmpdirname)
            meteor.tmp_dir = Path(tmpdirname)
            meteor.mapping_dir = Path(tmpdirname) / "map"
            meteor.fastq_dir = Path(tmpdirname)
            downloader = Downloader(meteor, Downloader.TEST_CATALOGUE, False, True)
            downloader.execute()
            fastq_importer = FastqImporter(meteor, meteor.tmp_dir, False, None)
            fastq_importer.execute()
            meteor.fastq_dir = Path(tmpdirname) / "test"
            meteor.ref_dir = meteor.ref_dir / "mock"
            counter = Counter(
                meteor, "total", "end-to-end", 80, 0.97, 100, 10, False, True
            )
            counter.execute()
            meteor.fastq_dir = Path(tmpdirname) / "test2"
            counter = Counter(
                meteor, "total", "end-to-end", 80, 0.97, 100, 10, False, True
            )
            counter.execute()
            meteor.mapped_sample_dir = meteor.mapping_dir / "test"
            meteor.strain_dir = Path(tmpdirname) / "strain"
            strain_detector = Strain(meteor, 100, 2, 2, 0.2, 1, 1, 0.2, 10, False)
            strain_detector.execute()
            meteor.mapped_sample_dir = meteor.mapping_dir / "test2"
            strain_detector = Strain(meteor, 100, 2, 2, 0.2, 1, 1, 0.2, 10, False)
            strain_detector.execute()
            meteor.tree_dir = Path(tmpdirname) / "tree"
            trees = TreeBuilder(
                meteor,
                0.1,
                4,
                800,
                600,
                None,
            )
            trees.execute()
    # Close logging
    logger.handlers[0].close()


if __name__ == "__main__":
    main()
