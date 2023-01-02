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

__author__ = "Amine Ghozlane (amine.ghozlane@inrae.fr)"
__version__  = "4.3"
__copyright__ = "GPLv3"
__date__ = "2022"

#-------------------------- MODULES IMPORTATION -------------------------------#
import sys
import os
import argparse
import configparser
import logging
from logging.handlers import RotatingFileHandler
import datetime
import glob
import re
import subprocess
# from dataclasses import dataclass

#---------------------------- CLASS DEFINITION --------------------------------#
class color:
   BOLD = '\033[1m'
   END = '\033[0m'

#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return os.path.abspath(path)


def isdir(path): # pragma: no cover
    """Check if path is an existing file.
      
    :param path: Path to the directory

    :raises ArgumentTypeError: If directory doesn't exist

    :return: (str) Path 
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return os.path.abspath(path) + os.sep


def get_log(path_log):
    """
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    # Create log file
    file_handler = RotatingFileHandler(path_log, 'a', 1000000, 1)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Stream in the the console
    ## TO REMOVE IF daemon
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)
    return logger


def get_arguments(): # pragma: no cover
    """
    Meteor help and arguments
    
    No arguments
    
    Return : parser
    """
    parser = argparse.ArgumentParser(description=color.BOLD + __doc__ + color.END)
    subparsers = parser.add_subparsers(title = 'positional arguments', help="Select activity")
    # # Mappping commands
    download_parser = subparsers.add_parser('download',
        help='Download catalog')
    reference_parser = subparsers.add_parser('build',
        help='Index reference')
    reference_parser.add_argument("-i", "--input", dest='input_file',
        type = str, required = True, help = "Input fasta filename.")
    reference_parser.add_argument("-p", "--'refRootDir'", dest='refRootDir',
        type = str, required = True, help = "Output path of the reference repository.")
    reference_parser.add_argument("-n", "--refName", dest='refName',
        type = str, required = True, help = "Name of the reference (ansi-string without space).")
    # reference_parser.add_argument('-1', '--no-bowtie1', dest="no_bowtie1", action='store_true', help='no index for bowtie1')
    # reference_parser.add_argument("-2", "--no-bowtie2", dest="no_bowtie2", action="store_true", help="no index for bowtie2")
    reference_parser.add_argument("-t", "--threads", dest='threads', default=4,
        type = int, help = "Thread count for bowtie2 (if available).")
    fastq_parser = subparsers.add_parser('fastq',
        help='Import fastq files')
    fastq_parser.add_argument("-i", "--fastqDir", dest='fastq_dir',
        type = isdir, required = True, help = """Path to a directory containing all input FASTQ files.
FASTQ files must have the extension .FASTQ or .FQ.
For paired-ends files must be named : 
    file_R1.[fastq/fq] & file_R2.[fastq/fq] 
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq]""")
    fastq_parser.add_argument("-p", "--projectName", dest='project_name',
        type = str, required = True, help = "Project name (ansi-string without space).")
    fastq_parser.add_argument("-m", "--maskSampleName", dest='mask_sample_name',
        type = str, required = True, help = "Regular expression for extracting sample name.")
    fastq_parser.add_argument("-t", "--techno", dest='techno', choices = ["illumina", "proton"],
        default="Illumina", required = True, type = str, help = "Sequencing technology (default: Illumina).")
    fastq_parser.add_argument("-c", "--compress", dest='iscompressed',
        default=False, action="store_true", help = "Fastq files are compressed.")
    fastq_parser.add_argument("-d", "--dispatch", dest='isdispatched',
        default=False, action="store_true", help = "Fastq files are already dispatched in directories.")
    mapping_parser = subparsers.add_parser('counter',
        help='Map reads against a gene catalog')
    mapping_parser.add_argument("-i", "--indir",
        type = isdir, required = True, action = "store", metavar = "DIR",
        help = """Path to input directory""")
    profiler_parser = subparsers.add_parser('profiler',
        help='Map reads against a gene catalog')
    return parser.parse_args(args=None if sys.argv[2:] else ['--help'])

#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Let's logging
    now = datetime.datetime.now()
    path_log = "meteor_" + now.strftime("%Y%m%d_%H%M") + ".log"
    # Get arguments
    args = get_arguments()
    logger = get_log(path_log)
    print(args)
    # Import FASTQ UGLY
    if hasattr(args, "fastq_dir"):
        if args.isdispatched:
            fastq_file_list = glob.glob(args.fastq_dir + "*" + os.sep + "*.f*q*")
        else:
            print(args.fastq_dir + "*.{fq,fastq}*")
            fastq_file_list = glob.glob('{}{}'.format(args.fastq_dir, "*.f*q*"))
        ext_R1 = ("1.fq", "1.fastq", "R1.fastq", "R1.fastq.gz", "R1.fq.gz")
        ext_R2 = ("2.fq", "2.fastq", "R2.fastq", "R2.fastq.gz", "R2.fq.gz")
        for fastq_file in fastq_file_list:
            print("Import ", fastq_file)
            logger.info("Import ", fastq_file)
            full_sample_name = os.path.basename(fastq_file)
            if args.iscompressed:
                full_sample_name = ".".join(full_sample_name.split(".")[:-1])
            print(full_sample_name)
            # Extract paired-end info
            tag = "single"
            if full_sample_name.endswith(ext_R1):
                tag = "1"
            elif full_sample_name.endswith(ext_R2):
                tag = "2"
            full_sample_name = ".".join(full_sample_name.split(".")[:-1])
            if args.isdispatched:
                sample_name = os.path.basename(os.path.dirname(fastq_file))
            else:
                # split full sample name (in fact library/run name) in order to extract sample_name according to regex mask
                full_sample_name_array = re.search(args.mask_sample_name, full_sample_name)
                # risk of TypeError: 'NoneType' object is not subscriptable
                sample_name = full_sample_name_array[0]
        
            if args.isdispatched:
                sample_dir = os.path.dirname(fastq_file) + os.sep
            else:
                # create directory for the sample and move fastq file into
                sample_dir = os.path.dirname(fastq_file) + os.sep + sample_name + os.sep
                if not os.path.isdir(sample_dir):
                    os.makedirs(sample_dir)
                os.rename(fastq_file, sample_dir + os.path.basename(fastq_file))
            config = configparser.ConfigParser()
            config["sample_info"] = {
                "sample_name": sample_name,
                "condition_name": "NA", # What is this ?
                "project_name" : args.project_name,
                "sequencing_date": "1900-01-01", # Then it is useless
                "sequencing_device": args.techno,
                "census_status": 0, # what is this ?
                "read_length": -1, # Then it is useless
                "tag": tag,
                "full_sample_name": full_sample_name
            }
            config["sample_file"] = {
                "fastq_file": os.path.basename(fastq_file),
                "is_compressed": args.iscompressed
            }
            with open(sample_dir + full_sample_name + "_census_stage_0.ini", 'w') as configfile:
                config.write(configfile)
    elif hasattr(args, "input_file"):
        # Create reference genome directory if it does not already exist
        ref_dir = os.path.join(args.refRootDir, args.refName)
        if not os.path.exists(ref_dir):
            os.makedirs(ref_dir)

        # Create subdirectories for fasta files and reference indices
        fasta_dir = os.path.join(ref_dir, "fasta")
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)

        database_dir = os.path.join(ref_dir, "database")
        if not os.path.exists(database_dir):
            os.makedirs(database_dir)
        # Read input fasta file and create new fasta file for each chromosome or contig
        output_annotation_file = os.path.join(database_dir, args.refName + '_lite_annotation')
        output_fasta_file = os.path.join(fasta_dir, args.refName + '.fasta')
        with open(args.input_file, 'rt') as input_fasta:
            with open(output_annotation_file, "wt") as output_annotation:
                with open(output_fasta_file, "wt") as output_fasta:
                    gene_id = 1
                    for line in input_fasta:
                        if line.startswith(">"):
                            output_annotation.write(line[1:].split(" ")[0].strip() +"\n")
                            output_fasta.write(">" + str(gene_id) + "\n")
                            gene_id += 1
                        else:
                            output_fasta.write(line)
        # Generate configuration file for reference genome
        config = configparser.ConfigParser()
        config["reference_info"] = {
            "reference_name": args.refName,
            "entry_type": "fragment", # Why ?
            "reference_date": datetime.datetime.now().strftime("%Y%m%d"),
            "database_type": "text",
            "HAS_LITE_INFO": 1
        }

        config["reference_file"] = {
            #IS_LARGE_REFERENCE_STR: 1,
            "database_dir": "database_dir", #WTF ?
            "fasta_dir": "fasta_dir", #WTF ?
            "fasta_file_count": 1,
            # is it possible to have several fasta
            "fasta_file_count": args.refName + '.fasta'
        }

        # if not args.no_bowtie2:
        index_prefix = os.path.join(database_dir, args.refName)
        subprocess.check_call(["bowtie2-build", '-f', '-t', str(args.threads), output_fasta_file, output_fasta_file] )
        config["bowtie2_index"] = {
            # "is_large_reference": "1", # WTF 
            "is_DNA_space_indexed": 1,
            "dna_space_bowtie_index_prefix_name_1": args.refName
        }

        # Write configuration file
        with open(os.path.join(ref_dir, args.refName + '_reference.ini'), 'wt') as config_file:
            config.write(config_file)


if __name__ == '__main__':
    main()