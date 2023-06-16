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

__version__ = "3.3"
__copyright__ = "GPLv3"
__date__ = "2023"

import sys
import logging
from argparse import ArgumentParser, ArgumentTypeError, Namespace
from pathlib import Path
from meteor.session import Component
from meteor.fastqimporter import FastqImporter
from meteor.referencebuilder import ReferenceBuilder
from meteor.counter import Counter
from meteor.downloader import Downloader
from meteor.profiler import Profiler
from tempfile import TemporaryDirectory


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
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise ArgumentTypeError(msg)
    return myfile


def isdir(path: str) -> Path:  # pragma: no cover
    """Check if path can be valid directory.

    :param path: Path to the directory

    :raises ArgumentTypeError: If directory does not exist

    :return: (str) Path object of the directory
    """
    mydir = Path(path)
    # if not mydir.is_dir():
    if mydir.is_file():
        msg = f"{mydir.name} is a file."
        raise ArgumentTypeError(msg)
        # else:
        #     msg = f"{mydir.name} does not exist."
    return mydir


def get_arguments() -> Namespace:  # pragma: no cover
    """
    Meteor help and arguments

    No arguments

    Return : parser
    """
    parser = ArgumentParser(description=Color.BOLD + __doc__ + Color.END, prog="Meteor")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(title="positional arguments",
                                       help="Select activity", dest="command",
                                       required=True)
    # Mappping commands
    download_parser = subparsers.add_parser("download", help="Download catalog")
    download_parser.add_argument("-i", dest="user_choice", type=str, required=True,
                                 choices=["chicken_caecal", "human_oral", "human_gut", "mouse_gut",
                                          "rabbit_gut", "rat_gut", "pig_gut"],
                                 help="Select the catalogue to download.")
    download_parser.add_argument("-c", dest="check_md5", action="store_true",
                                 help="Check the md5sum of the catalogue.")
    download_parser.add_argument("-fast", dest="taxonomy", action="store_true",
                                 help="Select the short catalogue version for taxonomical analysis.")
    download_parser.add_argument("-o", dest="ref_dir", type=isdir, required=True,
                                 help="Output directory.")
    reference_parser = subparsers.add_parser("build", help="Index reference")
    reference_parser.add_argument("-i", dest="input_fasta_file", type=isfile,
                                  required=True, help="Input fasta filename (compressed format accepted).")
    reference_parser.add_argument("-o", dest="ref_dir", type=isdir, required=True,
                                  help="Output path of the reference repository.")
    reference_parser.add_argument("-n", dest="ref_name", metavar="REFERENCE_NAME",
                                  type=str, required=True,
                                  help="Name of the reference (ansi-string without space).")
    reference_parser.add_argument("-t", dest="threads", default=1, type=int,
                                  help="Threads count.")
    reference_parser.add_argument("-no_pysam", dest="pysam_test", action="store_false",
                                  help="Execute original meteor")
    fastq_parser = subparsers.add_parser("fastq", help="Import fastq files")
    fastq_parser.add_argument("-i", dest="input_fastq_dir", type=isdir, required=True,
                              help="""Path to a directory containing all input fastq files.
                                        FASTQ files must have the extension .fastq or .fq.
                                        For paired-ends files must be named :
                                            file_R1.[fastq/fq] & file_R2.[fastq/fq]
                                                            or
                                            file_1.[fastq/fq] & file_2.[fastq/fq].
                                        If compressed, [gz,bz2,xz] are accepted.""")
    fastq_parser.add_argument("-p", dest="ispaired", default=False, action="store_true",
                              help="Fastq files are paired.")
    fastq_parser.add_argument("-m", dest="mask_sample_name", type=str,
                              help="Regular expression for extracting sample name.")
    fastq_parser.add_argument("-n", dest="project_name", type=str, required=True,
                              help="Project name (ansi-string without space).")
    fastq_parser.add_argument("-o", dest="fastq_dir", type=isdir, required=True,
                              help="Output path of the fastq repository.")
    # fastq_parser.add_argument("-c", dest="iscompressed", default=False,
    #     action="store_true", help = "Fastq files are compressed.")
    # fastq_parser.add_argument("-d", dest="isdispatched", default=False, action="store_true",
    #                           help="Fastq files are already dispatched in directories.")
    mapping_parser = subparsers.add_parser("mapping", help="Map reads against a gene catalog")
    mapping_parser.add_argument("-i", dest="fastq_dir", type=isdir, required=True,
                                help="""Path to sample directory, containing the sample sequencing metadata
                                        (files ending with _census_stage_0.ini)""")
    mapping_parser.add_argument("-r", dest="ref_dir", type=isdir, required=True,
                                help="Path to reference directory (Path containing *_reference.ini)")
    mapping_parser.add_argument("-o", dest="mapping_dir", type=isdir, required=True,
                                help="""Path to project directory, containing mapping
                                        and profile data (e.g. /projects/project_dir)""")
    mapping_parser.add_argument("-c", dest="counting_type", type=str,
                                default="smart_shared",
                                # "shared_reads",
                                choices=["total", "smart_shared", "unique", "best"],
                                help="Counting type string (default smart_shared_reads).")
    mapping_parser.add_argument("-p", dest="mapping_type", type=str,
                                choices=["local", "end-to-end"], default="end-to-end",
                                help="Counting type (Default end-to-end)")
    mapping_parser.add_argument("-trim", dest="trim", type=int, default=80,
                                help="Trim reads for mapping (default 80. If 0, no trim)")
    mapping_parser.add_argument("-id", dest="identity_threshold", type=float, default=0.95,
                                help="Aligned reads should have an identity to reference > 0.95 (default)."
                                "If 0, no filtering)")
    mapping_parser.add_argument("-align", dest="alignment_number", type=int, default=10000,
                                help="Number alignments considered for each read (default 10000)")
    mapping_parser.add_argument("-k", dest="keep_bam", action="store_true",
                                help="Save the bam files")
    mapping_parser.add_argument("-tmp", dest="tmp_path", type=isdir,
                                help="Path to the directory where temporary files (e.g. sam) are stored")
    mapping_parser.add_argument("-m", dest="mapping_only", action="store_true",
                                help="Execute mapping only")
    mapping_parser.add_argument("-n", dest="counting_only", action="store_true",
                                help="Execute counting only")
    mapping_parser.add_argument("-t", dest="threads", default=1, type=int,
                                help="Threads count.")
    mapping_parser.add_argument("-no_pysam", dest="pysam_test", action="store_false",
                                help="Execute original meteor")
    profiling_parser = subparsers.add_parser("profile", help="Performs profiling")
    profiling_parser.add_argument("-i", dest="mapping_dir", type=isdir, required=True,
                                help="""Path to project directory, containing mapping
                                        and profile data (e.g. /projects/project_dir)""")
    profiling_parser.add_argument("-r", dest="ref_dir", type=isdir, required=True,
                                help="Path to reference directory (Path containing *_reference.ini)")
    # profiling_parser.add_argument("-o", dest="profiling_dir", type=isdir, required=True,
    #                             help="""Path to project directory, containing mapping
    #                                     and profile data (e.g. /projects/project_dir)""")
    profiling_parser.add_argument("-n", dest="normalization", type=str,
                                  choices=["coverage", "rarefaction", "rpkm", "weighted_non_null_normalization"],
                                 help="Normalization techniques applied to gene abundances.")
    profiling_parser.add_argument("-l", dest="rarefaction_level", type=int,
                                  help="Rarefy according to a given level")
    profiling_parser.add_argument("-s", dest="statistics", type=str, choices=["specie", "ard", "kegg"], default="specie",
                                 help="""Group gene abundances into one given abundances (Default specie).""")
    profiling_parser.add_argument("-tmp", dest="tmp_path", type=isdir,
                                help="Path to the directory where temporary files (e.g. sam) are stored")
    subparsers.add_parser("test", help="Test meteor installation")
    return parser.parse_args(args=None if sys.argv[1:] else ["--help"])


def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Let us logging
    logger = get_logging()
    # Create a meteor dataset
    meteor = Component
    # Import FASTQ
    if args.command == "fastq":
        meteor.fastq_dir = args.fastq_dir
        # args.isdispatched,
        fastq_importer = FastqImporter(meteor, args.input_fastq_dir,
                                       args.ispaired,
                                       args.mask_sample_name,
                                       args.project_name)
        fastq_importer.execute()
    # Import reference
    elif args.command == "build":
        meteor.ref_name = args.ref_name
        meteor.ref_dir = args.ref_dir
        meteor.threads = args.threads
        reference_builder = ReferenceBuilder(meteor=meteor, input_fasta=args.input_fasta_file,
                                             pysam_test=args.pysam_test)
        reference_builder.execute()
    # Run mapping
    elif args.command == "mapping":
        meteor.fastq_dir = args.fastq_dir
        meteor.mapping_dir = args.mapping_dir
        meteor.ref_dir = args.ref_dir
        meteor.tmp_path = args.tmp_path
        meteor.threads = args.threads
        counter = Counter(meteor, args.counting_type, args.mapping_type,
                          args.trim, args. identity_threshold, args.alignment_number,
                          args.counting_only, args.mapping_only, args.keep_bam, args.pysam_test)
        counter.execute()
    # Run download catalogues
    elif args.command == "download":
        meteor.ref_name = args.user_choice
        meteor.ref_dir = args.ref_dir
        if args.taxonomy:
            args.user_choice += "_taxo"
        downloader = Downloader(meteor, args.user_choice, args.check_md5)
        downloader.execute()
    elif args.command == "profile":
        meteor.mapping_dir = args.mapping_dir
        meteor.ref_dir = args.ref_dir
        meteor.tmp_path = args.tmp_path
        profiler = Profiler(meteor, args.normalization, args.rarefaction_level, args.statistics)
        profiler.execute()
    # Testing
    else:
        with TemporaryDirectory() as tmpdirname:
            meteor.ref_name = "test"
            meteor.ref_dir = Path(tmpdirname) / "ref"
            meteor.threads = 1
            meteor.tmp_path = Path("")
            meteor.tmp_dir = Path(tmpdirname)
            meteor.mapping_dir = Path(tmpdirname) / "map"
            meteor.fastq_dir = Path(tmpdirname)
            downloader = Downloader(meteor, "test", True)
            downloader.execute()
            fastq_importer = FastqImporter(meteor, meteor.tmp_dir, False, None, "test_project")
            fastq_importer.execute()
            meteor.fastq_dir = Path(tmpdirname) / "test"
            counter = Counter(meteor, "best", "end-to-end", 80, 1, False, False, False, True)
            counter.execute()
    # Close logging
    logger.handlers[0].close()


if __name__ == "__main__":
    main()
