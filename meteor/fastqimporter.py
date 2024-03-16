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

"""Import and prepare fastq"""

import sys
import logging
import re
from itertools import product, chain
from pathlib import Path
from dataclasses import dataclass, field
from typing import Iterator
from meteor.session import Session, Component


@dataclass
class FastqImporter(Session):
    """FastqImporter handle the fastq import"""

    meteor: type[Component]
    input_fastq_dir: Path
    ispaired: bool
    mask_sample_name: str | None
    ext_r1: tuple = field(default_factory=tuple)
    ext_r2: tuple = field(default_factory=tuple)
    ext: tuple = field(default_factory=tuple)

    def __post_init__(self) -> None:
        self.ext_r1 = tuple(self.extension("1"))
        self.ext_r2 = tuple(self.extension("2"))
        self.ext = tuple(self.short_extension())

    def extension(self, pair: str) -> Iterator[str]:  # pragma: no cover
        """Get all possible extensions for a given fastq including pairing info

        :param pair: A string giving the strand
        :return: (Generator) A generator of a string from a given combination
        """
        for ext in product(
            self.meteor.sequence, pair, self.meteor.extension, self.meteor.compression
        ):
            yield "".join(ext)

    def short_extension(self) -> Iterator[str]:  # pragma: no cover
        """Get all possible extensions for a given fastq

        :return: (Generator) A generator of a string from a given combination
        """
        for short_ext in product(self.meteor.extension, self.meteor.compression):
            yield "".join(short_ext)

    def replace_ext(self, fastq_filename: str) -> str:
        """Replace all fastq/compressed extension to get a fullpathname

        :param fastq_filename: Name of the fastq file
        :return: (str) A string without the expected extension in the name
        """
        for e in self.ext:
            fastq_filename = fastq_filename.replace(e, "")
        return fastq_filename

    def get_paired_dirname(self, fastq_filename: str, tag: str) -> str:
        """Replace all fastq/compressed extension and pairing to get a sample_name

        :param fastq_filename: Name of the fastq file
        :return: (str) A string without the expected extension in the name
        """
        for e in tuple(self.extension(tag)):
            fastq_filename = fastq_filename.replace(e, "")
        return fastq_filename

    def get_fastq_files(self) -> Iterator[Path]:  # pragma: no cover
        """Find all fastq files in the given input"""
        files_expected_ext = chain.from_iterable(
            self.input_fastq_dir.glob("*" + e) for e in self.short_extension()
        )
        return (f for f in files_expected_ext if f.is_file())

    def get_tag(self, fastq_filename: str) -> str | None:
        """Extract paired-end info

        :param fastq_filename: Name of the fastq file
        :return: (str) A string giving pairing status
        """
        if fastq_filename.endswith(self.ext_r1):
            return "1"
        elif fastq_filename.endswith(self.ext_r2):
            return "2"
        else:
            return None

    def set_fastq_config(
        self, sample_name: str, tag: str, fastq_file: Path, full_sample_name: str
    ) -> dict:  # pragma: no cover
        """Set configuration for fastq

        :param sample_name: Sample name
        :param tag: Identified tag (single or 1 or 2)
        :param fastq_file: A fastq file path
        :param full_sample_name: A fastq file name without the extension
        :return: (Dict) A dict configuration
        """
        config = {
            "meteor_version": self.meteor.version,
            "sample_info": {
                "sample_name": sample_name,
                "tag": tag,
                "full_sample_name": full_sample_name,
            },
            "sample_file": {"fastq_file": fastq_file.name},
        }
        return config

    def execute(self) -> None:
        """Dispatch the fastq file"""
        logging.info("Starting import of fastq files from %s", self.input_fastq_dir)
        fastq_files = list(self.get_fastq_files())
        if not fastq_files:
            logging.error("No fastq file detected")
            sys.exit(1)
        else:
            logging.info("%d fastq files detected", len(fastq_files))

        num_imported_fastq_files = 0
        samples_names = set()
        for fastq_file in fastq_files:
            # Get rid of all possible extension
            full_sample_name = self.replace_ext(fastq_file.name)
            if self.ispaired:
                # Extract paired-end info
                tag = self.get_tag(fastq_file.name)
                if not tag:
                    logging.error(
                        "Pairing tag (1 or 2) is not detected in %s", fastq_file
                    )
                    sys.exit(1)
            else:
                tag = "single"
            # split full sample name (in fact library/run name) in order
            if self.mask_sample_name:
                # to extract sample_name according to regex mask
                full_sample_name_array = re.search(
                    self.mask_sample_name, full_sample_name
                )
                if full_sample_name_array:
                    sample_name = full_sample_name_array[0]
                else:
                    logging.warning("Regular expression does not match %s", fastq_file)
                    continue
            else:
                if self.ispaired:
                    sample_name = self.get_paired_dirname(fastq_file.name, tag)
                else:
                    sample_name = full_sample_name
            logging.info("Importing %s in sample %s", fastq_file, sample_name)
            # Create directory for the sample and symlink fastq file into
            samples_names.add(sample_name)
            sample_dir = self.meteor.fastq_dir / sample_name
            sample_dir.mkdir(exist_ok=True, parents=True)
            sym_fastq = Path(sample_dir / fastq_file.name)
            if not sym_fastq.is_symlink():
                sym_fastq.symlink_to(fastq_file.resolve())
            # Create a configuration
            config_fastq = self.set_fastq_config(
                sample_name, tag, sym_fastq, full_sample_name
            )
            config_path = sample_dir / f"{full_sample_name}_census_stage_0.json"
            self.save_config(config_fastq, config_path)
            num_imported_fastq_files += 1
        logging.info(
            "%d/%d fastq files imported in %d samples",
            num_imported_fastq_files,
            len(fastq_files),
            len(samples_names),
        )
