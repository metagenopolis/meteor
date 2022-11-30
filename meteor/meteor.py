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
from argparse import RawTextHelpFormatter

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
    return path


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
    return path


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
    fastq_parser = subparsers.add_parser('fastq',
        help='Import fastq files, Useless ?')
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
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()