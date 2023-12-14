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

"""Find msp files, align and """
from pathlib import Path
import argparse
from collections import defaultdict


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
        raise argparse.ArgumentTypeError(msg)
        # else:
        #     msg = f"{mydir.name} does not exist."
    return mydir


def main():
    parser = argparse.ArgumentParser(
        description="Example script with an integer parameter in a specific range"
    )

    # Define an integer parameter with a range between 1 and 100
    parser.add_argument(
        "-i",
        dest="input_dir",
        type=isdir,
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "-o",
        dest="output_dir",
        type=isdir,
        required=True,
        help="Output directory.",
    )

    args = parser.parse_args()
    msp_file = defaultdict(list)
    for filepath in args.input_dir.glob("**/" + "msp_*"):
        msp_file[filepath.name].append(filepath)

    # Concatenate files that occur in more than one directory
    for filename, paths in msp_file.items():
        if len(paths) > 1:
            res = args.output_dir / f"{filename}_concatenated"
            with res.open("w") as outfile:
                for path in paths:
                    with open(path, "r") as infile:
                        outfile.write(infile.read())


if __name__ == "__main__":
    main()
