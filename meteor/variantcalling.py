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

"""Effective variant calling"""

import logging
import sys
import gzip
import lzma
import bgzip
import pickle
import shutil
import pandas as pd
import numpy as np
from subprocess import CalledProcessError, run, Popen, PIPE
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from meteor.session import Session, Component
from time import perf_counter
from tempfile import NamedTemporaryFile
from packaging.version import parse
from pysam import (
    AlignmentFile,
    AlignmentHeader,
    FastaFile,
    VariantFile,
    faidx,
    index,
    tabix_index,
    bcftools,
)
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from typing import ClassVar
from tqdm import tqdm


def depth_array(
    cram: AlignmentFile,
    gene_name: str,
    gene_length: int,
    fasta: FastaFile,
    max_depth: int,
) -> np.ndarray:
    """Per-position read count of a gene, deletions and ref-skips excluded.

    Same pileup and counting rule as VariantCalling.count_reads_in_gene, but
    returned as a dense array of length gene_length.
    """
    depth = np.zeros(gene_length, dtype=np.int64)
    for pileupcolumn in cram.pileup(
        contig=gene_name,
        start=0,
        end=gene_length,
        stepper="all",
        max_depth=max_depth,
        fastafile=fasta,
        multiple_iterators=False,
    ):
        pos = pileupcolumn.reference_pos
        if pos < gene_length:
            depth[pos] = sum(
                True
                for pileupread in pileupcolumn.pileups
                if not pileupread.is_del and not pileupread.is_refskip
            )
    return depth


def runs_below(
    depth: np.ndarray, min_depth: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Maximal runs of equal depth over [0, len) whose depth is < min_depth.

    Vectorised equivalent of VariantCalling.group_consecutive_positions.
    """
    length = len(depth)
    change = np.flatnonzero(depth[1:] != depth[:-1]) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [length]))
    values = depth[starts]
    keep = values < min_depth
    return starts[keep], ends[keep], values[keep]


def low_cov_worker(
    cram_file: str,
    reference_file: str,
    genes: list[tuple[str, int]],
    max_depth: int,
    min_depth: int,
) -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """Low coverage runs for a batch of genes (runs in a worker process)."""
    result = []
    with AlignmentFile(cram_file, "rc", reference_filename=reference_file) as cram:
        with FastaFile(filename=reference_file) as fasta:
            for gene_name, gene_length in genes:
                starts, ends, values = runs_below(
                    depth_array(cram, gene_name, gene_length, fasta, max_depth),
                    min_depth,
                )
                if len(starts) > 0:
                    result.append((gene_name, starts, ends, values))
    return result


def low_cov_to_dict(low_cov_sites: "pd.DataFrame | dict") -> dict:
    """{gene index value: int array (n, 2) of [start, end) intervals}.

    Replaces per-gene `.loc` lookups on a non-unique, unsorted index, which
    scan the whole table each time (O(rows) per gene).
    """
    if isinstance(low_cov_sites, dict):
        return low_cov_sites
    if len(low_cov_sites) == 0:
        return {}
    codes, uniques = pd.factorize(low_cov_sites.index)
    order = np.argsort(codes, kind="stable")
    intervals = low_cov_sites[["startpos", "endpos"]].to_numpy()[order]
    bounds = np.cumsum(np.bincount(codes, minlength=len(uniques)))[:-1]
    return dict(zip(uniques, np.split(intervals, bounds)))


def gene_ignore_to_dict(gene_ignore: "pd.DataFrame | dict") -> dict:
    """{gene_id: gene_length} of genes fully replaced by gaps."""
    if isinstance(gene_ignore, dict):
        return gene_ignore
    if len(gene_ignore) == 0 or "gene_length" not in gene_ignore.columns:
        return {}
    return dict(zip(gene_ignore.index, gene_ignore["gene_length"].astype(int)))


def run_freebayes_chunk(
    temp_ref_file_path: str,
    bed_chunk_file: Path,
    cram_file: Path,
    vcf_chunk_file: Path,
    min_snp_depth: int,
    min_frequency: float,
    ploidy: int,
    tmp_dir: Path,
):
    """Function to run freebayes on a chunk of the BED file (i.e., a portion of the genome)."""
    try:
        with Popen(
            [
                "freebayes",
                "-i",  # no indel
                "-u",  # no complex observation that may include ins
                "--pooled-continuous",
                "--haplotype-length",
                str(0),
                "--min-alternate-count",
                str(1),
                "--min-coverage",
                str(min_snp_depth),
                "--min-alternate-fraction",
                str(min_frequency),
                "--min-mapping-quality",
                str(0),
                "--use-duplicate-reads",
                "-t",
                str(
                    bed_chunk_file.resolve()
                ),  # BED region chunk for parallel execution
                "-p",
                str(ploidy),
                "-f",
                temp_ref_file_path,  # Path to temporary reference FASTA file
                "-b",
                str(cram_file.resolve()),  # BAM/CRAM alignment file
            ],
            stdout=PIPE,
            stderr=PIPE,
        ) as freebayes_process:
            # Capture output from freebayes
            freebayes_output, freebayes_error = freebayes_process.communicate()
            if freebayes_error:
                logging.error(
                    "Error processing chunk %s: %s", bed_chunk_file, freebayes_error
                )
                return None
            elif freebayes_process.returncode == 0:
                # Compress output using bgzip
                with vcf_chunk_file.open("wb") as raw:
                    with bgzip.BGZipWriter(raw) as fh:
                        fh.write(freebayes_output)
                bcftools.sort('-Oz', '-o', str(vcf_chunk_file.resolve()), '-T', str(tmp_dir), str(vcf_chunk_file.resolve()), catch_stdout=False)
                tabix_index(str(vcf_chunk_file.resolve()), preset="vcf", force=True)
            else:
                logging.error(
                    "Freebayes process failed for chunk %s, return code: %d",
                    bed_chunk_file,
                    freebayes_process.returncode,
                )
                return None
        return vcf_chunk_file

    except CalledProcessError as e:
        logging.error(
            "Freebayes failed for chunk %s with return code %d",
            bed_chunk_file,
            e.returncode,
        )
        logging.error("Error output: %s", e.output)
        raise e


@dataclass
class VariantCalling(Session):
    """Run freebayes"""

    # from https://www.bioinformatics.org/sms/iupac.html
    IUPAC: ClassVar[dict] = {
        ("A",): "A",
        ("T",): "T",
        ("G",): "G",
        ("C",): "C",
        ("A", "G"): "R",
        ("C", "T"): "Y",
        ("C", "G"): "S",
        ("A", "T"): "W",
        ("G", "T"): "K",
        ("A", "C"): "M",
        ("C", "G", "T"): "B",
        ("A", "G", "T"): "D",
        ("A", "C", "T"): "H",
        ("A", "C", "G"): "V",
        ("A", "C", "G", "T"): "N",
    }

    meteor: type[Component]
    census: dict
    max_depth: int
    min_depth: int
    min_snp_depth: int
    min_frequency: float
    ploidy: int
    core_size: int

    # freebayes BED chunks per thread (dynamic load balancing)
    FREEBAYES_CHUNKS_PER_THREAD: ClassVar[int] = 4

    def gene_counts(self) -> dict[int, float]:
        """Raw gene counts (count table of the mapping step)"""
        if getattr(self, "_gene_counts", None) is None:
            counts = pd.read_csv(
                self.matrix_file,
                sep="\t",
                names=["gene_id", "gene_length", "value"],
                header=0,
                compression="xz",
            )
            counts = counts[counts["value"] > 0]
            self._gene_counts = dict(zip(counts["gene_id"], counts["value"]))
        return self._gene_counts

    def set_variantcalling_config(
        self,
        cram_file: Path,
        vcf_file: Path,
        consensus_file: Path,
        freebayes_version: str,
    ) -> dict:  # pragma: no cover
        """Define the census 1 configuration

        :param cmd: A string of the specific parameters
        :param cram_file: A path to the sam file
        :return: (Dict) A dict object with the census 1 config
        """
        config = {
            "meteor_version": self.meteor.version,
            "sample_info": self.census["census"]["sample_info"],
            "mapping": {
                "reference_name": self.census["reference"]["reference_info"][
                    "reference_name"
                ],
                "cram_name": cram_file.name,
            },
            "variant_calling": {
                "variant_calling_tool": "freebayes",
                "variant_calling_version": freebayes_version,
                "variant_calling_date": datetime.now().strftime("%Y-%m-%d"),
                "vcf_name": vcf_file.name,
                "consensus_name": consensus_file.name,
                "min_snp_depth": str(self.min_snp_depth),
                "min_frequency": str(self.min_frequency),
                "max_depth": str(self.max_depth),
            },
        }
        return config

    def create_bed_chunks(
        self, merged_df: pd.DataFrame, num_chunks: int, tmp_dir: Path
    ) -> list[Path]:
        """
        Divide the merged_df DataFrame into `num_chunks` BED chunks,
        with each chunk containing multiple `msp_name` groups.
        Each chunk will be written to a separate temporary BED file.
        """
        bed_chunks = []

        # Get unique msp_name groups
        msp_names = merged_df["msp_name"].unique()
        total_msp_names = len(msp_names)

        # Calculate how many MSPs each chunk should include
        base_chunk_size = total_msp_names // num_chunks

        # Number of chunks that will get an extra MSP name (due to remainder)
        remainder = total_msp_names % num_chunks

        start_idx = 0  # Starting index for slicing the msp_names

        # Split the msp_names into balanced chunks
        for i in range(num_chunks):
            # Determine the chunk size: add 1 to the base size if within the remainder limit
            chunk_size = base_chunk_size + 1 if i < remainder else base_chunk_size

            # Determine the end index for this chunk
            end_idx = start_idx + chunk_size

            # Select the subset of MSPs for this chunk
            current_msp_names = msp_names[start_idx:end_idx]

            # Subset the DataFrame for only these selected MSPs
            chunk_df = merged_df[merged_df["msp_name"].isin(current_msp_names)]

            # Create a temporary BED file for this chunk
            temp_bed_file = NamedTemporaryFile(suffix=".bed", dir=tmp_dir, delete=False)

            # Write the chunk DataFrame subset to the temporary BED file
            chunk_df[["gene_id", "startpos", "gene_length"]].to_csv(
                temp_bed_file.name, sep="\t", index=False, header=False
            )

            # Add the path to the temporary file to the list of bed_chunks
            bed_chunks.append(Path(temp_bed_file.name))

            # Update the starting index for the next chunk
            start_idx = end_idx

        return bed_chunks  # Return the list of file paths.

    def create_balanced_bed_chunks(
        self,
        merged_df: pd.DataFrame,
        gene_weight: dict[int, float],
        num_chunks: int,
        tmp_dir: Path,
    ) -> list[Path]:
        """Split MSPs into `num_chunks` BED files of similar expected cost.

        The cost of an MSP is the number of reads counted on its core genes.
        MSPs are assigned greedily (longest-processing-time first) and the
        chunks are returned from the most to the least loaded, so that the
        process pool starts with the slowest chunks.
        """
        weights = (
            merged_df.assign(
                weight=merged_df["gene_id"].map(gene_weight).fillna(0.0) + 1.0
            )
            .groupby("msp_name", sort=False)["weight"]
            .sum()
            .sort_values(ascending=False, kind="stable")
        )
        num_chunks = max(1, min(num_chunks, len(weights)))
        loads = [0.0] * num_chunks
        members: list[list[str]] = [[] for _ in range(num_chunks)]
        for msp_name, weight in weights.items():
            i = min(range(num_chunks), key=loads.__getitem__)
            loads[i] += weight
            members[i].append(msp_name)
        bed_chunks = []
        for i in sorted(range(num_chunks), key=lambda k: -loads[k]):
            chunk_df = merged_df[merged_df["msp_name"].isin(set(members[i]))]
            temp_bed_file = NamedTemporaryFile(suffix=".bed", dir=tmp_dir, delete=False)
            chunk_df[["gene_id", "startpos", "gene_length"]].to_csv(
                temp_bed_file.name, sep="\t", index=False, header=False
            )
            bed_chunks.append(Path(temp_bed_file.name))
        return bed_chunks

    def write_marker_reference(
        self, reference_file: Path, gene_ids: list[int], tmp_dir: Path
    ) -> str:
        """Uncompressed FASTA (+ .fai) restricted to the genes given to freebayes.

        The filtered CRAM only holds alignments on these genes, so freebayes
        never needs the rest of the catalogue.
        """
        with NamedTemporaryFile(
            suffix=".fasta", dir=tmp_dir, delete=False, mode="wt"
        ) as out:
            with FastaFile(filename=str(reference_file.resolve())) as fasta:
                for gene_id in gene_ids:
                    out.write(f">{gene_id}\n{fasta.fetch(str(gene_id))}\n")
            path = out.name
        faidx(path)
        return path

    def write_marker_alignments(
        self,
        cram_file: Path,
        reference_file: Path,
        gene_ids: set[int],
        tmp_dir: Path,
    ) -> tuple[Path, list[int]]:
        """Indexed BAM of the filtered alignments whose header only lists `gene_ids`
        (plus any gene carrying an alignment).

        The filtered CRAM inherits the bowtie2 header (one @SQ per catalogue
        gene, ~10M for hs_10_4_gut). htslib/freebayes materialise that header
        in every process: ~16 GB RSS and ~90 s per freebayes call on
        hs_10_4_gut, against ~0.2 GB and ~1 s with the marker-only header
        (identical VCF records).

        :return: (BAM path, gene ids of the new header in header order)
        """
        with AlignmentFile(
            str(cram_file.resolve()),
            "rc",
            reference_filename=str(reference_file.resolve()),
            threads=self.meteor.threads,
        ) as cram:
            names = cram.references
            lengths = cram.lengths
            present = self.reference_ids_in_cram(cram_file)
            keep = [
                tid
                for tid, name in enumerate(names)
                if tid in present or int(name) in gene_ids
            ]
            new_tid = np.full(len(names), -1, dtype=np.int64)
            new_tid[keep] = np.arange(len(keep))
            header = AlignmentHeader.from_dict(
                {
                    "HD": {"VN": "1.6", "SO": "coordinate"},
                    "SQ": [{"SN": names[t], "LN": lengths[t]} for t in keep],
                }
            )
            bam_path = Path(
                NamedTemporaryFile(suffix=".bam", dir=tmp_dir, delete=False).name
            )
            with AlignmentFile(
                str(bam_path), "wb", header=header, threads=self.meteor.threads
            ) as out:
                for read in cram:
                    # coordinate order is kept: tids are renumbered monotonically
                    read.reference_id = int(new_tid[read.reference_id])
                    if read.next_reference_id >= 0:
                        read.next_reference_id = int(new_tid[read.next_reference_id])
                    out.write(read)
        index(str(bam_path))
        return bam_path, [int(names[t]) for t in keep]

    def reference_ids_in_cram(self, cram_file: Path) -> set[int]:
        """Reference ids having alignments, read from the CRAM index (.crai)"""
        crai = Path(f"{cram_file}.crai")
        if not crai.exists():
            index(str(cram_file.resolve()))
        ref_ids: set[int] = set()
        with gzip.open(crai, "rt") as index_fh:
            for line in index_fh:
                ref_id = int(line.split("\t", 1)[0])
                if ref_id >= 0:
                    ref_ids.add(ref_id)
        return ref_ids

    def group_consecutive_positions(
        self, position_count_dict: dict, gene_name: str, gene_length: int
    ):
        """Runs of equal read count below min_depth, as a DataFrame"""
        depth = np.zeros(gene_length, dtype=np.int64)
        for pos, count in position_count_dict.items():
            if 0 <= pos < gene_length:
                depth[pos] = count
        starts, ends, values = runs_below(depth, self.min_depth)
        return pd.DataFrame(
            {
                "gene_id": [gene_name] * len(starts),
                "startpos": starts,
                "endpos": ends,
                "coverage": values,
            },
            columns=["gene_id", "startpos", "endpos", "coverage"],
        )

    def count_reads_in_gene(
        self,
        cram: AlignmentFile,
        gene_name: str,
        gene_length: int,
        Fasta: FastaFile,
    ):
        """
        Counts the number of reads at each position in a specified gene from a BAM file,
        ignoring reads with gaps or deletions at that position.

        Parameters:
        cram (pysam.AlignmentFile): A pysam AlignmentFile object representing the BAM file
        gene_name (str): The name of the gene for which reads are to be counted

        Returns:
        dict: A dictionary with positions (int) as keys and the number of reads (int) at each position as values.
        """
        reads_dict: dict[int, int] = defaultdict(int)
        # gene_length = bam.lengths[bam.references.index(gene_name)]
        # Extraction of positions having at least one read
        for pileupcolumn in cram.pileup(
            contig=gene_name,
            start=0,
            end=gene_length,
            stepper="all",
            # stepper="samtools",
            max_depth=self.max_depth,
            fastafile=Fasta,
            multiple_iterators=False,
        ):
            read_count = sum(
                True
                for pileupread in pileupcolumn.pileups
                if not pileupread.is_del and not pileupread.is_refskip
            )
            reads_dict[pileupcolumn.reference_pos] = read_count

        return reads_dict

    def merge_vcf_files(self, vcf_file_list, output_vcf) -> None:
        """Merge variant records (handling the same positions in multiple VCFs)."""
        variant_dict = defaultdict(list)

        # Collect all records from all files
        for i, vcf_file in enumerate(vcf_file_list):
            with VariantFile(vcf_file, threads=self.meteor.threads) as vcf_in:
                #  Get the header from the first input VCF
                if i == 0:
                    vcf_header = vcf_in.header
                for rec in vcf_in:
                    # Use (chrom, pos) tuple as key to merge records based on positions
                    variant_dict[(rec.chrom, rec.pos)].append(rec)
        # Write the merged VCF output
        with VariantFile(
            str(output_vcf.resolve()),
            "w",
            header=vcf_header,
            threads=self.meteor.threads,
        ) as vcf_out:
            for _, rec_list in variant_dict.items():
                vcf_out.write(rec_list[0])

    # @memory_profiler.profile
    def filter_low_cov_sites(
        self,
        cram_file: Path,
        reference_file: Path,
        gene_subset: set[int] | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Report, per gene, the runs of positions below the depth threshold

        :param cram_file:   Path to the input cram file
        :param reference_file: Path to the catalogue fasta
        :param gene_subset: Only these genes are analysed (the genes written in
            the consensus). The filtered CRAM has no alignment elsewhere, so
            other genes would be piled up for nothing.
        :return: (low coverage runs indexed by gene_id as str,
                  genes below min_depth indexed by gene_id)
        """
        # Get list of genes
        gene_interest = pd.read_csv(
            self.matrix_file,
            sep="\t",
            names=[
                "gene_id",
                "gene_length",
                "coverage",
            ],
            header=0,
            compression="xz",
        )
        # Round the coverage column to 0 decimal places
        gene_interest["coverage"] = gene_interest["coverage"].round(0)

        # Convert the 'coverage' column to a sparse column using SparseDtype
        gene_interest["coverage"] = gene_interest["coverage"].astype(
            pd.SparseDtype(int, 0)
        )
        any_covered_gene = bool((gene_interest["coverage"] >= self.min_depth).any())
        if gene_subset is not None:
            gene_interest = gene_interest[gene_interest["gene_id"].isin(gene_subset)]
        # We need to do more work for these genes
        # their total count is above the min_depth
        # but we need to find on which positions they do.
        gene_tofilter = gene_interest[gene_interest["coverage"] >= self.min_depth][
            ["gene_id", "gene_length"]
        ]
        genes = [
            (str(gene_id), int(gene_length))
            for gene_id, gene_length in gene_tofilter.itertuples(index=False)
        ]
        cram_path = str(cram_file.resolve())
        ref_path = str(reference_file.resolve())
        workers = max(1, self.meteor.threads)
        if workers == 1 or len(genes) < 2 * workers:
            results = [
                low_cov_worker(cram_path, ref_path, genes, self.max_depth, self.min_depth)
            ]
        else:
            # small batches keep the pool balanced; order of results is kept
            batch = max(1, min(2000, len(genes) // (workers * 8)))
            batches = [genes[i : i + batch] for i in range(0, len(genes), batch)]
            with ProcessPoolExecutor(max_workers=workers) as executor:
                results = list(
                    executor.map(
                        low_cov_worker,
                        [cram_path] * len(batches),
                        [ref_path] * len(batches),
                        batches,
                        [self.max_depth] * len(batches),
                        [self.min_depth] * len(batches),
                    )
                )
        runs = [run for result in results for run in result]
        # All these genes are going to be replaced by gaps
        # Their count is below the threshold level
        gene_ignore = gene_interest[
            gene_interest["coverage"] < self.min_depth
        ].set_index("gene_id")

        if len(runs) == 0 and (gene_subset is None or not any_covered_gene):
            logging.error("No low coverage regions detected, it might be linked to no coverage at all")
            sys.exit(1)
        if len(runs) == 0:
            sum_cov_bed = pd.DataFrame(
                columns=["gene_id", "startpos", "endpos", "coverage"]
            ).set_index("gene_id")
        else:
            sum_cov_bed = pd.DataFrame(
                {
                    "gene_id": np.repeat(
                        np.array([run[0] for run in runs], dtype=object),
                        [len(run[1]) for run in runs],
                    ),
                    "startpos": np.concatenate([run[1] for run in runs]),
                    "endpos": np.concatenate([run[2] for run in runs]),
                    "coverage": np.concatenate([run[3] for run in runs]),
                }
            ).set_index("gene_id")
        return sum_cov_bed, gene_ignore

    # @memory_profiler.profile
    def create_consensus(
        self,
        reference_file,
        consensus_file,
        low_cov_sites,
        gene_ignore,
        vcf_file,
        bed_file,
    ):
        """Generate a consensus sequence by applying VCF variants to the provided reference genome."""
        # Read the CSV file using pandas
        bed_set = sorted(
            set(
                pd.read_csv(bed_file, usecols=[0], sep="\t", header=None)
                .iloc[:, 0]
                .astype(int)
            )
        )
        # low_cov_sites / gene_ignore may be DataFrames (pickled by older runs,
        # tests) or dicts: index them once instead of a `.loc` scan per gene
        low_cov = low_cov_to_dict(low_cov_sites)
        ignore_length = gene_ignore_to_dict(gene_ignore)
        gap = self.meteor.DEFAULT_GAP_CHAR
        gap_byte = gap.encode()
        with VariantFile(str(vcf_file.resolve()), threads=self.meteor.threads) as vcf:
            with FastaFile(filename=str(reference_file.resolve())) as Fasta:
                with lzma.open(consensus_file, "wt", preset=0) as consensus_f:
                    # Iterate over all reference sequences in the fasta file
                    for gene_id in tqdm(
                        bed_set, desc="Creating consensus", unit="gene"
                    ):
                        ref = str(gene_id)
                        if gene_id in ignore_length:
                            consensus_f.write(f">{gene_id}\n")
                            consensus_f.write(gap * int(ignore_length[gene_id]) + "\n")
                            continue
                        consensus = np.frombuffer(
                            Fasta.fetch(ref).encode("ascii"), dtype="S1"
                        ).copy()
                        # Apply variants from VCF
                        for record in vcf.fetch(ref):
                            ##INFO=<ID=RO,Number=1,Type=Integer,Description="Count of full observations of the reference haplotype.">
                            ##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
                            reference_frequency = record.info["RO"] / (
                                record.info["RO"] + np.sum(record.info["AO"])
                            )
                            if reference_frequency >= self.min_frequency:
                                keep_alts = tuple(sorted(list(record.alleles)))
                            else:
                                keep_alts = tuple(sorted(list(record.alts)))
                            max_len = max(map(len, keep_alts))
                            # MNV vase
                            if max_len > 1:
                                for i in range(max_len):
                                    mnv = tuple(
                                        sorted(
                                            set(
                                                keep_alts[k][i]
                                                for k in range(len(keep_alts))
                                            )
                                        )
                                    )
                                    consensus[record.start + i] = self.IUPAC[mnv]
                            else:
                                consensus[record.start] = self.IUPAC[keep_alts]
                        # Mark low coverage positions as uncertain
                        intervals = low_cov.get(ref)
                        if intervals is not None:
                            for startpos, endpos in intervals:
                                consensus[startpos:endpos] = gap_byte
                        consensus_f.write(f">{gene_id}\n")
                        consensus_f.write(consensus.tobytes().decode("ascii") + "\n")
                        del consensus

    def execute(self) -> None:
        """Call variants reads"""
        # Start mapping
        cram_file = (
            self.census["mapped_sample_dir"]
            / f"{self.census['census']['sample_info']['sample_name']}.cram"
        )
        self.matrix_file = (
            self.census["mapped_sample_dir"]
            / f"{self.census['census']['sample_info']['sample_name']}.tsv.xz"
        )
        vcf_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.vcf.gz"
        )
        low_cov_sites_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}.pickle"
        )
        consensus_file = (
            self.census["directory"]
            / f"{self.census['census']['sample_info']['sample_name']}_consensus.fasta.xz"
        )
        reference_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["fasta_dir"]
            / self.census["reference"]["reference_file"]["fasta_filename"]
        )
        msp_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["msp"]["filename"]
        )
        # print(self.census)
        annotation_file = (
            self.meteor.ref_dir
            / self.census["reference"]["reference_file"]["database_dir"]
            / self.census["reference"]["annotation"]["gene_id"]["filename"]
        )
        msp_content = self.load_data(msp_file)
        gene_details = self.load_data(annotation_file)
        freebayes_exec = run(
            ["freebayes", "--version"],
            check=False,
            capture_output=True,
        )
        if freebayes_exec.returncode != 0:
            logging.error(
                "Checking freebayes failed:\n%s", freebayes_exec.stderr.decode("utf-8")
            )
            sys.exit(1)
        freebayes_version = (
            freebayes_exec.stdout.decode("utf-8").strip().split("  ")[1][1:]
        )
        if parse(freebayes_version) < self.meteor.MIN_FREEBAYES_VERSION:
            logging.error(
                "The freebayes version %s is outdated for meteor. Please update freebayes to >= %s.",
                freebayes_version,
                self.meteor.MIN_FREEBAYES_VERSION,
            )
            sys.exit(1)

        start = perf_counter()
        startfreebayes = perf_counter()
        temp_ref_file_path = None
        # Prepare the gene data by merging content and creating the necessary fields for the BED format
        msp_content = msp_content[msp_content["gene_category"] == "core"]
        msp_content = (
            msp_content.groupby("msp_name")
            .head(self.core_size)
            .reset_index(drop=True)
        )  # Limit to core_size per `msp_name`

        # Merge with gene details
        merged_df = pd.merge(msp_content, gene_details, on="gene_id")

        # Add BED columns (we assume `startpos` is 0 and `gene_length` is the length of the gene)
        merged_df["startpos"] = 0
        merged_df["gene_length"] = merged_df["gene_length"].astype(
            int
        )  # Ensure these are integers
        result_df = merged_df[["gene_id", "startpos", "gene_length"]]
        temp_bed_file = NamedTemporaryFile(suffix=".bed", dir=self.meteor.tmp_dir, delete=False)
        result_df.to_csv(temp_bed_file, sep="\t", index=False, header=False)
        bed_genes = set(merged_df["gene_id"].astype(int))
        marker_bam: Path | None = None
        if not (vcf_file.exists() and low_cov_sites_file.exists()):
            # Alignment file + FASTA restricted to marker genes: the filtered CRAM
            # header lists the whole catalogue, which every freebayes / pysam
            # process would otherwise load (no full catalogue decompression either)
            startmarker = perf_counter()
            marker_bam, marker_genes = self.write_marker_alignments(
                cram_file, reference_file, bed_genes, self.meteor.tmp_dir
            )
            temp_ref_file_path = self.write_marker_reference(
                reference_file, marker_genes, self.meteor.tmp_dir
            )
            logging.info(
                "Marker alignments/reference (%d genes) written in %f seconds",
                len(marker_genes),
                perf_counter() - startmarker,
            )
        if vcf_file.exists():
            logging.info("Vcf already exist, skipping freebayes..")
            vcf_chunk_files = []
        else:
            logging.info("Run freebayes")
            assert marker_bam is not None and temp_ref_file_path is not None
            # Create bed_chunk files. Each file stores multiple `msp_name` regions,
            # several chunks per thread, balanced on the reads counted per MSP
            bed_chunks = self.create_balanced_bed_chunks(
                merged_df,
                self.gene_counts(),
                self.meteor.threads * self.FREEBAYES_CHUNKS_PER_THREAD,
                self.meteor.tmp_dir,
            )
            # List to store the VCF chunk files
            vcf_chunk_files = [
                NamedTemporaryFile(
                    suffix=".vcf.gz", dir=self.meteor.tmp_dir, delete=False
                ).name
                for _ in bed_chunks
            ]
            # Use ProcessPoolExecutor to run freebayes in parallel on each BED chunk
            failed_chunks = []
            with ProcessPoolExecutor(max_workers=self.meteor.threads) as executor:
                futures = {
                    executor.submit(
                        run_freebayes_chunk,
                        temp_ref_file_path,  # Pass the path to the reference file
                        bed_chunk_file,  # Each BED chunk
                        marker_bam,
                        Path(vcf_chunk_file),
                        self.min_snp_depth,
                        self.min_frequency,
                        self.ploidy,
                        self.meteor.tmp_dir,
                    ): bed_chunk_file
                    for bed_chunk_file, vcf_chunk_file in zip(
                        bed_chunks, vcf_chunk_files
                    )
                }

                # Iterate through completed futures
                for future in as_completed(futures):
                    bed_chunk = futures[future]
                    try:
                        vcf_chunk_file = future.result()
                        logging.info(
                            "Processed BED chunk %s -> VCF chunk %s",
                            bed_chunk,
                            vcf_chunk_file,
                        )
                        if vcf_chunk_file is None:
                            failed_chunks.append(bed_chunk)
                    except Exception as exc:
                        logging.error("Error processing chunk %s: %s", bed_chunk, exc)
                        failed_chunks.append(bed_chunk)
            if failed_chunks:
                logging.error(
                    "freebayes failed on %d/%d chunks, variant calling aborted",
                    len(failed_chunks),
                    len(bed_chunks),
                )
                sys.exit(1)

            logging.info("All chunks have been processed")
            # Combine VCF chunk files into the final VCF
            if len(vcf_chunk_files) > 1:
                logging.info("Merging vcf")
                self.merge_vcf_files(vcf_chunk_files, vcf_file)
            else:
                shutil.copyfile(str(vcf_chunk_files[0]), str(vcf_file.resolve()))
        logging.info(
            "Completed freebayes step in %f seconds", perf_counter() - startfreebayes
        )
        # Index the vcf file
        startindexing = perf_counter()
        if not Path(f"{vcf_file}.tbi").exists():
            logging.info("Indexing")
            bcftools.sort('-Oz','-o', str(vcf_file.resolve()), '-T', str(self.meteor.tmp_dir), str(vcf_file.resolve()), catch_stdout=False)
            tabix_index(str(vcf_file.resolve()), preset="vcf", force=True)
        else:
            logging.info("Index already exist, skipping...")
        logging.info(
            "Completed indexing step in %f seconds", perf_counter() - startindexing
        )
        # The columns of the tab-delimited BED file are also CHROM, POS
        # and END (trailing columns are ignored), but coordinates are
        # 0-based, half-open. To indicate that a file be treated as BED
        # rather than the 1-based tab-delimited file, the file must have
        # the ".bed" or ".bed.gz" suffix (case-insensitive).
        startlowcovpython = perf_counter()
        if low_cov_sites_file.exists():
            logging.info("Loading low coverage regions")
            with low_cov_sites_file.open("rb") as file:
                # Load the data from the file
                data = pickle.load(file)
            low_cov_sites = data["low_cov_sites"]
            gene_ignore = data["gene_ignore"]
        else:
            logging.info("Detecting low coverage regions")
            assert marker_bam is not None and temp_ref_file_path is not None
            low_cov_sites, gene_ignore = self.filter_low_cov_sites(
                marker_bam,
                Path(temp_ref_file_path),
                gene_subset=bed_genes,
            )
            # Open a file for writing the pickle data (binary write mode)
            with low_cov_sites_file.open("wb") as file:
                # Dump the data into the file
                pickle.dump(
                    {"low_cov_sites": low_cov_sites, "gene_ignore": gene_ignore}, file
                )
        logging.info(
            "Completed low coverage regions filtering step in %f seconds",
            perf_counter() - startlowcovpython,
        )
        logging.info("Consensus creation")
        startconsensuspython = perf_counter()
        self.create_consensus(
            reference_file,
            consensus_file,
            low_cov_sites,
            gene_ignore,
            vcf_file,
            temp_bed_file.name,
        )
        logging.info(
            "Completed consensus step in %f seconds",
            perf_counter() - startconsensuspython,
        )
        logging.info("Completed SNP calling in %f seconds", perf_counter() - start)
        config = self.set_variantcalling_config(
            cram_file, vcf_file, consensus_file, freebayes_version
        )
        self.save_config(config, self.census["Stage3FileName"])
        # Cleanup temporary files
        temporary_files = (
            [temp_ref_file_path] + vcf_chunk_files + [f"{temp_ref_file_path}.fai"]
            if temp_ref_file_path is not None
            else vcf_chunk_files
        )
        if marker_bam is not None:
            temporary_files += [str(marker_bam), f"{marker_bam}.bai"]
        for temp_file in temporary_files:
            p = Path(temp_file)
            if p.exists():
                p.unlink(missing_ok=True)
