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

__version__  = "3.3"
__copyright__ = "GPLv3"
__date__ = "2022"

import sys
import logging
from argparse import ArgumentParser, ArgumentTypeError, Namespace
# import urllib.request
from pathlib import Path
from MeteorFastqImporter import FastqImporter
from MeteorReferenceBuilder import ReferenceBuilder
from MeteorCounter import Counter

class Color:
    BOLD = "\033[1m"
    END = "\033[0m"


def get_logging():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    print("Here please")
    # Create log file
    file_handler = logging.FileHandler("meteor.log", 'a', encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Stream in the the console
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)
    return logger


def isfile(path: str)-> Path: # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file doesn't exist

    :return: (Path) Path
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise ArgumentTypeError(msg)
    return myfile


def isdir(path: str) -> Path: # pragma: no cover
    """Check if path is an existing file.

    :param path: Path to the directory

    :raises ArgumentTypeError: If directory doesn't exist

    :return: (str) Path
    """
    mydir = Path(path)
    # if not mydir.is_dir():
    if mydir.is_file():
        msg = f"{mydir.name} is a file."
        raise ArgumentTypeError(msg)
        # else:
        #     msg = f"{mydir.name} does not exist."
    return mydir



# def download_url(url:str, output_path: str):
#     with DownloadProgressBar(unit='B', unit_scale=True,
#                              miniters=1, desc=url.split('/')[-1]) as t:
#         urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def get_arguments()->Namespace: # pragma: no cover
    """
    Meteor help and arguments

    No arguments

    Return : parser
    """
    parser = ArgumentParser(description=Color.BOLD + __doc__ + Color.END)
    parser.add_argument("-t", dest='threads', default=4,
        type = int, help = "Threads count.")
    subparsers = parser.add_subparsers(title = 'positional arguments', help="Select activity", dest= "command")
    # # Mappping commands
    download_parser = subparsers.add_parser('download',
        help='Download catalog')
    reference_parser = subparsers.add_parser('build',
        help='Index reference')
    reference_parser.add_argument("-i", dest='input_fasta_file',
        type = isfile, required = True, help = "Input fasta filename.")
    reference_parser.add_argument("-o", dest='reference_dir',
        type = isdir, required = True, help = "Output path of the reference repository.")
    reference_parser.add_argument("-n", dest='ref_name',
        type = str, required = True, help = "Name of the reference (ansi-string without space).")
    fastq_parser = subparsers.add_parser('fastq', help='Import fastq files')
    fastq_parser.add_argument("-i", dest='fastq_dir', type = isdir, 
        help = """Path to a directory containing all input fastq files.
FASTQ files must have the extension .fastq or .fq.
For paired-ends files must be named :
    file_R1.[fastq/fq] & file_R2.[fastq/fq]
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq].
If compressed, [gz,bz2,xz] are accepted.""", required = True)
    fastq_parser.add_argument("-p", dest='project_name',
        type = str, required = True, help = "Project name (ansi-string without space).")
    fastq_parser.add_argument("-m", dest='mask_sample_name',
        type = str, required = True, help = "Regular expression for extracting sample name.")
    # fastq_parser.add_argument("-t", dest='techno', choices = ["illumina", "proton"],
    #     default="Illumina", type = str, help = "Sequencing technology (default: Illumina).")
    fastq_parser.add_argument("-c", dest='iscompressed',
        default=False, action="store_true", help = "Fastq files are compressed.")
    fastq_parser.add_argument("-d", dest='isdispatched',
        default=False, action="store_true", help = "Fastq files are already dispatched in directories.")
    mapping_parser = subparsers.add_parser('mapping',
        help='Map reads against a gene catalog')
    # mapping_parser.add_argument("-w", dest="workflow_ini", type=isfile, required=True,
        # help="Path to meteor configuration file, e.g. workflow.ini")
    mapping_parser.add_argument("-i", dest="input_dir", type = isdir, required = True,
        help = """Path to sample directory, containing the sample sequencing metadata
        (files ending with _census_stage_0.ini)""")
    mapping_parser.add_argument("-r", dest="reference_dir", type=isdir, required=True,
        help="""Path to sample directory, containing the sample sequencing metadata
        (files ending with _census_stage_0.ini)""")
    mapping_parser.add_argument("-o", dest="mapping_dir", type=isdir, required=True,
        help="path to project directory, containing mapping and profile data (e.g. /projects/project_dir)")
    mapping_parser.add_argument("-c", dest="counting_type", type=str, nargs='+', default=["smart_shared_reads"],
        choices = ["total_reads", "shared_reads", "smart_shared_reads", "unique_reads"],
    help="""counting type string (default smart_shared_reads.""")
    # mapping_parser.add_argument("-s", dest="strand", type=str, choices=["direct", "undirect"],
    #     help="(where direct/undirect means direct/undirect strand)")
    # mapping_parser.add_argument("-a", dest="coverage", type=str, choices=["coverage", "avg_coverage"],
    #     help="coverage or average coverage")
    mapping_parser.add_argument("-t", dest="tmp_dir", type=isdir,
        help="path to the directory where temporary files (e.g. sam) are stored")
    # mapping_parser.add_argument("-e", dest="clean_tmp_path", action="store_true",
    #     help="Remove the temporary files directory after the task counting is finished")
    # mapping_parser.add_argument("-l", dest="remove_lock", action="store_true",
    #     help="ignore and remove any lock file")
    mapping_parser.add_argument("-f", dest="force", action="store_true",
        help="force overwriting, ignoring previous results and lock files")
    mapping_parser.add_argument("-m", dest="mapping_only", action="store_true",
        help="execute mapping only")
    mapping_parser.add_argument("-n", dest="counting_only", action="store_true",
        help="execute counting only")
    return parser.parse_args(args=None if sys.argv[2:] else ['--help'])


#==============================================================
# Main program
#==============================================================
def main()->None: # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Let us logging
    logger = get_logging()
    # logging_meteor = LoggingMeteorSession(MeteorSession)
    # Import FASTQ
    if args.command == "fastq":
        fastq_importer = FastqImporter(args.threads,
            args.isdispatched, args.fastq_dir,
            args.mask_sample_name, args.project_name)
        fastq_importer.execute()
    # Import reference
    elif args.command ==  "build":
        reference_builder = ReferenceBuilder(args.threads, 
            args.input_fasta_file, args.reference_dir, args.ref_name)
        reference_builder.execute()
    elif args.command == "mapping":
        counter = Counter(args.threads, args.input_dir,
        args.mapping_dir, args.counting_type, args.reference_dir, args.tmp_dir,
        args.force, args.counting_only, args.mapping_only)
        counter.execute()
    #     m= MeteorSession(logger, args.tmp_path)
    #     # counting_type = args.counting_type
    #     m.ProcessJob(args.workflow_ini,
    #         args.project_path, args.input_path,
    #         args.force,
    #         args.mapping_dir, args.counting_type,
    #         args.counting_only, args.mapping_only)
    # Close logging
    logger.handlers[0].close()


if __name__ == '__main__':
    main()