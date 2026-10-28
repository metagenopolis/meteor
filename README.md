# Meteor

[![install with conda](https://img.shields.io/conda/vn/bioconda/meteor?color=green&label=bioconda%2Fmeteor&logo=anaconda)](https://anaconda.org/bioconda/meteor)
![Anaconda downloads](https://anaconda.org/bioconda/meteor/badges/downloads.svg)
[![PyPI](https://img.shields.io/pypi/v/METEOR?label=pypi%20package)](https://pypi.org/project/meteor/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/meteor?period=total&units=INTERNATIONAL_SYSTEM&left_color=GRAY&right_color=BLUE&left_text=meteor%20on%20PyPi)](https://pepy.tech/projects/meteor)
[![dockerhub](https://img.shields.io/docker/v/aghozlane/meteor?label=aghozlane/meteor&logo=docker)](https://hub.docker.com/r/aghozlane/meteor/)
![Github Actions](https://github.com/metagenopolis/meteor/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/metagenopolis/meteor/graph/badge.svg?token=AXAEIUY7DX)](https://codecov.io/gh/metagenopolis/meteor)


## Introduction

Meteor is a plateform for quantitative metagenomics profiling of complex ecosystems.
Meteor relies on genes catalogue to perform species-level taxonomic profiling (Bacteria, Archaea and Eukaryotes), functional analysis and strain-level population structure inference.


## Dependencies

Besides python packages dependencies, Meteor requires:
- python 3.11+
- [bowtie2 2.3.5+](https://github.com/BenLangmead/bowtie2)
- [freebayes 1.3.6+](https://github.com/freebayes/freebayes)

## Installation

Meteor is available with conda which includes all its dependencies:
```
conda create --name meteor -c conda-forge -c bioconda meteor
```

Or with pip with a recent Python 3.10+:
```
pip install meteor
```
You can test the installation of meteor with:
```
meteor test
```

## Nextflow wrapper

For automated pipeline execution, a Nextflow wrapper `nf-meteor.nf` is available that streamlines the entire Meteor workflow:

```bash
nf-meteor.nf --in <fastq_dir> --catalogue_name <catalogue_name> --out <output_dir> --cpus <nb_cpus> -w <temp_work_dir>
Parameters:
--in Directory containing paired fastq.gz files (default: ).
--out Output directory (default: ).
--cpus Number of cpus to use (default: 4).
--catalogue_name Name of the prebuilt catalogue to use (default: none). Allowed values are: fc_1_3_gut, gg_13_6_caecal, clf_1_0_gut, hs_10_4_gut, hs_8_4_oral, hs_2_9_skin, mm_5_0_gut, oc_5_7_gut, rn_5_9_gut, ssc_9_3_gut
--catalogue Path to a catalogue (overrides --catalogue_name if both are provided).
--fast Enable fast mode for meteor (no functional analysis) (default: false).
--check_catalogue Check md5sum of the catalogue is compatible with the input reads (default: false).
```

## Getting started

A basic usage of meteor will require to:
1. **Download or build a reference catalogue**
2. **Structure the raw fastq files**
3. **Map reads against the reference catalogue**
4. **Compute taxonomical and/or functional abundances**
5. **Strain profiling**

### 1. Download a reference
-------------------------------------------

Meteor requires to download locally a microbial gene catalogue specif, either in 'full' or 'light' version. The 'full' version contains all genes of the catalogue, whereas the 'light' version contains only the marker genes that will be used to infer species abundance profiles. Of note, no functional profiling can be performed when using the 'light' version of a catalogue.

Ten catalogues are currently available:

|  Microbial gene catalogue | \<name\> | Genes count (M) | Metagenomic Species Pan-genomes (MSPs) |Size (full) (GB) | Size (light) (GB)  | Description  |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
|  *Canis lupus familiaris* | clf_1_0_gut  | 0.95  | 234 | 1.4 | 0.1 |[link](https://zenodo.org/records/16982647)
|  *Felis catus* | fc_1_3_gut  | 1.3  | 344 | 2.0 | 0.2 |[link](https://zenodo.org/records/16982752)
|  *Gallus gallus domesticus* | gg_13_6_caecal  | 13.6  | 2420 | 19.6 | 1.2 |[link](https://zenodo.org/records/16982807)
| *Homo sapiens gut* |  hs_10_4_gut | 10.4  | 1990 | 12.6 | 0.7 |[link](https://zenodo.org/records/16982835)
| *Homo sapiens oral*  |  hs_8_4_oral | 8.4  | 853 | 13.7 | 0.5 |[link](https://zenodo.org/records/16983006)
| *Homo sapiens skin*  |  hs_2_9_skin | 2.9  | 392 | 3.9 | 0.2 |[link](https://zenodo.org/records/16982908)
| *Mus musculus*  | mm_5_0_gut  | 5.0  | 1252 | 10.3 | 0.6 |[link](https://zenodo.org/records/16983064)
| *Oryctolagus cuniculus* | oc_5_7_gut  | 5.7 | 1053 | 8.0 | 0.4 |[link](https://zenodo.org/records/16983124)
| *Rattus norvegicus* | rn_5_9_gut  | 5.9 | 1627 | 7.0 | 0.6 |[link](https://zenodo.org/records/16983154)
| *Sus domesticus* | ssc_9_3_gut  | 9.3  | 1523 | 11.3 | 0.7 |[link](https://zenodo.org/records/16983194)

These references can be downloaded with the following command:
```
meteor download -i <name> -c -o <refdir>
```
The 'light' catalogues are available with the tag (--fast) :
```
meteor download -i <name> -c --fast -o <refdir>
```

### 2. Import fastq
-------------------
Meteor requires a first of fastq indexing:
```
meteor fastq -i <fastqdir>  [-p paired reads] -o <outputdir>
```
When multiple sequencing are available for a library, the option -m allows to group these samples.
Example:

Illumina_lib1-**SAMPLE_01**.fastq <br />
Illumina_lib1-**SAMPLE_02**.fastq <br />
Illumina_lib2-**SAMPLE_01**.fastq <br />
Illumina_lib2-**SAMPLE_02**.fastq <br />

In this case, the following command will group these samples the same library:
```
meteor fastq -i ./  -m SAMPLE_\\d+ -o outputdir
```

### 3. Mapping
----------------
The fastq files are mapped against a catalogue to generate a gene count table with the following command:
```
meteor mapping -i <fastqdir/sampledir> -r <refdir> -o <mappingdir>
```
We recommend to first filter out reads with low-quality, length < 60nt or belonging to the host.

### 4. Taxonomic and functional profiling
-------------------------

Genes from the catalogue are clustered into Metagenomic Species Pangeomes (MSP) with [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712), and are functionnaly annotated against [KEGG r107](https://academic.oup.com/nar/article/36/suppl_1/D480/2507484), [DBcan](https://academic.oup.com/nar/article/51/W1/W115/7147496?login=true) (carbohydate active enzymes) and [MUSTARD](https://www.nature.com/articles/s41564-018-0292-6) (antibiotic resistant determinants).

 MSP and functional profiles are computed from the gene count table with the following command:

```
meteor profile -i <mappingdir/sampledir> -o <profiledir> -r <refdir> -n coverage
```

The "-n" parameter ensures read count normalization for gene length. If omitted, no normalization will be performed on the gene table.

This profiling step will generate:
- a Species abundance table;
- an ARD abundance table (full catalogue only);
- a DBCAN abundance table (full catalogue only);
- a Gut Metabolic Modules ([GMM](https://www.nature.com/articles/nmicrobiol201688)) abundance table (from the KO annotation) (full catalogue only).
- a Gut Brain Modules ([GBM](https://www.nature.com/articles/s41564-018-0337-x)) abundance table (from the KO, EGGNOG and TIGRFAM annotations) (full catalogue only).

### 5. Merging

To merge output from different samples into a single table, use the following command:

```
meteor merge -i <profiledir> -r <refdir> -o <mergingdir>
```

All tables are merged by default, except the gene tables (`-g`). `-t` merges only the tables given, e.g. the species table alone to look for cross-sample contamination with CroCoDeEL:

```
meteor merge -i <profiledir> -r <refdir> -o <mergingdir> -t msp
```

Tables: `msp` (species), `genes` (gene counts, same as `-g`), `modules_completeness`, `modules`, and the annotation tables `mustard`, `resfinder`, `resfinderfg`, `dbcan`, `kegg` (each `_as_genes_sum` and `_as_msp_sum`). Functional tables need a complete catalogue.

`-l <samples.txt>` merges only the samples named in a text file, one sample name per line (first column if the file is tab-separated; empty lines and lines starting with `#` are ignored), in the order of the file. Samples of the list that are not found are reported in the log.

With `-m`, the species abundance table is also saved in the style of a merged MetaPhlAn table (`<prefix>_msp_mpa.tsv`), for the tools that read MetaPhlAn profiles. Each row is a clade from superkingdom (`k__`) to MSP (`t__`), named by its lineage (`k__Bacteria|p__Pseudomonadota|...|s__Escherichia_coli|t__msp_0005`). Values are relative abundances in percent of the sample: the MSPs written in `<prefix>_msp.tsv` (after the `-a`, `-n` and `-s` filters) sum to 100, and a clade is the sum of its MSPs, so each rank sums to 100. Ranks missing from the catalogue taxonomy are named `unclassified`. `--tax_lev` keeps a single rank (`k`, `p`, `c`, `o`, `f`, `g`, `s` or `t`). The first line names the meteor version, the catalogue and the GTDB release of its taxonomy.

### 5. Strain profiling
-------------------------

Meteor is capable of profiling strains in large metagenomic datasets. It identifies specific mutations from strains and applies them to the  gene catalog MSPs.

To use Meteor for strain profiling, use the following command:
```
meteor strain -i <mappingdir/sampledir> -o <straindir> -r <refdir>
```

Meteor computes mutation rates and trees between strains from samples using a GTR+GAMMA model with the following command:
```
meteor tree -i <straindir> -o <treedir>
```

`-l <samples.txt>` builds the trees with only the samples of the list (same format as for `meteor merge`). When the output directory already holds the tree of an MSP computed from other samples, it is recomputed.

This profiling step will generate:
- a mutation rate matrix;
- a fasta file of each strain for each sample;
- a table giving detailed comparison per strain. The file is a tab-separated values (TSV) file with one row per sample pair. Each row contains the following columns:

| Column | Description |
|--------|-------------|
| `sample1` | First sample in the comparison |
| `sample2` | Second sample in the comparison |
| `total_length` | Total length of the alignment (in bases) |
| `overlap_noN_info_count` | Number of positions where **both** samples have minimal information (A,C,G,T or IUPAC codes *excluding* N, gaps, and ?) |
| `overlap_noIUPAC_info_count` | Number of positions where **both** samples have maximal information (strictly A,C,G,T only) |
| `overlap_noN_info_pc` | Percentage of `total_length` with minimal information overlap |
| `overlap_noIUPAC_info_pc` | Percentage of `total_length` with maximal information overlap |
| `noN_info_pc_sample1` | In sample1, percentage of positions with minimal information |
| `noN_info_pc_sample2` | In sample2, percentage of positions with minimal information |
| `noIUPAC_info_pc_sample1` | In sample1, percentage of positions with maximal information |
| `noIUPAC_info_pc_sample2` | In sample2, percentage of positions with maximal information |
| `distance` | Genetic distance between samples (0.0 = identical) |
| `distance_category` | Categorical classification: `same_strain`, `same_subspecies`, or `divergent` |


Samples are automatically classified based on their genetic distance:

| Category | Distance Threshold | Approximate Similarity | Biological Interpretation |
|----------|-------------------|------------------------|---------------------------|
| `same_strain` | ≤ 0.0001 | ≥ 99.99% | Same strain/clone |
| `same_subspecies` | ≤ 0.015 | ≥ 97% | Same subspecies |
| `divergent` | > 0.015 | < 97% | Different lineage |

**Note**: 
- **Maximal Information**: Unambiguous nucleotides (A, C, G, T only)
- **Minimal Information**: Includes IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V) but excludes N, gaps (-), and unknown (?)

Samples sitting on abnormally long branches of a strain tree (rogue taxa) inflate their distances to every other sample. `meteor tree --treeshrink` runs the per-gene test of [TreeShrink](https://github.com/uym2/TreeShrink) (Mai & Mirarab, BMC Genomics 2018) on each MSP tree and writes, next to the tree:
- `<msp>_qcflags.tsv`: for each sample of the best removal sets, its diameter ratio (largest reduction of the tree diameter obtained by removing it), the threshold of the MSP (log-lognormal, alpha = 0.05) and whether it is flagged;
- `<msp>_shrink.tsv`: the distance matrix without the flagged samples (only when samples are flagged).

#### Strain sharing

At the end of `meteor tree`, the comparison tables of all MSPs are gathered into three tables:

- `strain_sharing_events.tsv`: one row per strain-sharing event, i.e. a sample pair sharing the strain of an MSP (`msp_name`, `sample1`, `sample2`, `relation`, `distance`, `overlap_noN_info_pc`, `threshold`);
- `strain_sharing_samples.tsv`: one row per sample pair, with the number of MSPs compared (`msp_compared`), shared (`msp_shared`) and the strain-sharing rate (`sharing_rate` = shared / compared);
- `strain_sharing_msp.tsv`: one row per MSP, with its threshold and, over all pairs, the related pairs and the unrelated pairs, the number of pairs compared and shared and their rate (`sharing_rate`, `transmissibility`, `unrelated_sharing_rate`).

Two samples are compared on an MSP when the sites informative in both cover at least `--min_overlap` % of the alignment (default 10). They share its strain when their distance is at most the threshold of the MSP.

With `--groups <groups.tsv>`, a tab-separated file giving the group of each sample (household, mother-infant pair, subject followed over time...), samples of the same group are related and samples of different groups unrelated:
- the threshold of each MSP is the 3rd percentile of the distances between unrelated samples (at least 50 of them), as the strain identity thresholds of StrainPhlAn (Valles-Colomer et al., Nature 2023); otherwise it is 0.0001 (`threshold_source` = `fixed`);
- `transmissibility` is the sharing rate of the related pairs, to compare with `unrelated_sharing_rate`, the background of the MSP.

Without `--groups`, every threshold is 0.0001 and no transmissibility is computed. Rates of an MSP are given from `--min_pairs` compared pairs (default 10).

A tree directory computed in GTR by meteor 2.0.22 or earlier has wrong distances in its comparison tables (issue #133): `meteor tree` recomputes them when it runs again in that directory.

The tree and the distance matrix are not modified. The test is implemented in meteor (no TreeShrink installation needed) and flags the same samples as `run_treeshrink.py -m per-gene` (TreeShrink 1.4.0). Existing trees are not rebuilt, so `meteor tree --treeshrink` on an existing tree directory only adds the QC files.

## Options of each command

The options of every command are also listed by `meteor <command> -h`.

### meteor download

| Option | Description | Default |
|---|---|---|
| `-i USER_CHOICE` | Select the catalogue to download. Choices: `fc_1_3_gut`, `gg_13_6_caecal`, `clf_1_0_gut`, `hs_10_4_gut`, `hs_8_4_oral`, `hs_2_9_skin`, `mm_5_0_gut`, `oc_5_7_gut`, `rn_5_9_gut`, `ssc_9_3_gut`. | required |
| `--fast` | Select the short catalogue version (only for taxonomic profiling). | off |
| `-o REF_DIR` | Directory where the downloaded catalogue is saved. | required |
| `-c` | Check the md5sum of the catalogue after download. | off |

### meteor fastq

| Option | Description | Default |
|---|---|---|
| `-i INPUT_FASTQ_DIR` | Directory containing all input fastq files with .fastq or .fq. extensions (gzip, bzip2 and xz compression accepted). Paired-end files must be named : file_R1.[fastq/fq] & file_R2.[fastq/fq] or file_1.[fastq/fq] & file_2.[fastq/fq] | required |
| `-p` | Fastq files are paired. | off |
| `-m MASK_SAMPLE_NAME` | Regular expression (between quotes) for extracting sample name. | none |
| `-o FASTQ_DIR` | Directory where the fastq repository is created. | required |

### meteor mapping

| Option | Description | Default |
|---|---|---|
| `-i FASTQ_DIR` | Directory corresponding to the sample to process. (contains sequencing metadata files ending with _census_stage_0.json) | required |
| `-r REF_DIR` | Directory corresponding to the gene catalog against which reads are mapped. (contains a file ending with *_reference.json) | required |
| `-o MAPPING_DIR` | Directory where mapping and raw gene counts of the sample are saved. | required |
| `-p MAPPING_TYPE` | Strategy to map reads against the catalogue. Choices: `end-to-end`, `local`. | `end-to-end` |
| `--trim TRIM` | Trim reads exceeding TRIM bases before mapping. If 0, no trim. | `80` |
| `--align ALIGNMENT_NUMBER` | Maximum number alignments to report for each read. | `100` |
| `-c COUNTING_TYPE` | Strategy to calculate raw gene counts. Choices: `total`, `smart_shared`, `unique`. | `smart_shared` |
| `--core_size CORE_SIZE` | Number of signature genes per species (MSP) used to estimate their respective abundance. | `100` |
| `--id IDENTITY_THRESHOLD` | Keep only the read alignments with a nucleotide identity >= IDENTITY_THRESHOLD. If 0.0, no filtering. | auto: `0.95` (complete catalogue), `0.98` (taxonomic catalogue) |
| `--ka` | Keep raw bowtie2 output in cram format. Required for calculating gene counts with another strategy. | off |
| `--kf` | Keep filtered alignments on marker genes. Required for strain analysis. | off |
| `--tmp TMP_PATH` | Directory where temporary files (e.g. cram) are stored | none |
| `-t THREADS` | Number of alignment threads to launch. | `1` |

### meteor profile

| Option | Description | Default |
|---|---|---|
| `-i MAPPED_SAMPLE_DIR` | Directory with raw gene counts of the sample to process. (contains a metadata file ending with _census_stage_1.json) | required |
| `-r REF_DIR` | Directory corresponding to the catalog used to generate raw gene counts. (contains a file ending with *_reference.json) | required |
| `-o PROFILED_SAMPLE_DIR` | Directory where species and functional abundance tables of the sample are saved. | required |
| `-l RAREFACTION_LEVEL` | Rarefaction level. If <= 0, no rarefation is performed. | `0` |
| `--seed SEED` | Seed of the random number generator used for rarefaction. | `1234` |
| `-n NORMALIZATION` | Normalization applied to raw gene counts. Choices: `coverage`, `fpkm`, `raw`. | `coverage` |
| `-c COVERAGE_FACTOR` | Multiplication factor for coverage normalization. | `100.0` |
| `--core_size CORE_SIZE` | Number of signature genes per species (MSP) used to estimate their respective abundance. | `100` |
| `--msp_filter MSP_FILTER` | Minimal proportion of core genes detected in a sample to consider a species (MSP) as present. | auto: `0.1` (complete catalogue), `0.2` (taxonomic catalogue) |
| `--completeness COMPLETENESS` | Cutoff above which a module is considered as present in a species. Value between 0.0 and 1.0. | `0.9` |

### meteor merge

| Option | Description | Default |
|---|---|---|
| `-i PROFILE_DIR` | Directory containing subdirectories (one per sample) with abundance tables to be merged. (each subdirectory contains a metadata file ending with _census_stage_2.json) | required |
| `-r REF_DIR` | Directory corresponding to the gene catalog used to generate the abundance tables. (contains a file ending with *_reference.json) | required |
| `-a MIN_MSP_ABUNDANCE` | Minimum msp abundance. | `0.0` |
| `-n MIN_MSP_OCCURRENCE` | Report only species (MSPs) occuring in at least n samples. | `1` |
| `-s` | Remove samples with no detected species (MSPs). | off |
| `-m` | Save the merged species abundance table in the style of MetaPhlAn: relative abundances (%) of the clades from superkingdom (k__) to MSP (t__). | off |
| `-b` | Save the merged species abundance table in biom format. | off |
| `--tax_lev MPA_TAXONOMIC_LEVEL` | Taxonomic level of the MetaPhlAn-style table (`-m`): `a` all levels, `k` kingdoms, `p` phyla, `c` classes, `o` orders, `f` families, `g` genera, `s` species, `t` MSPs only. | `a` |
| `-o MERGING_DIR` | Directory where the merged abundance tables are saved. | required |
| `-p PREFIX` | Prefix added to output filenames. | `output` |
| `-g` | Merge gene abundance tables. | off |
| `-t TABLE [TABLE ...]` | Merge only these tables: `msp` (species), `genes` (gene counts, same as `-g`), `modules_completeness`, `modules`, and the annotation tables `mustard`, `resfinder`, `resfinderfg`, `dbcan`, `kegg` (each `_as_genes_sum` and `_as_msp_sum`). Functional tables need a complete catalogue. | all but `genes` |
| `-l SAMPLE_FILE` | Merge only the samples of this file, in its order: one sample name per line (first column if tab-separated; empty lines and lines starting with # ignored). | none |

### meteor strain

| Option | Description | Default |
|---|---|---|
| `-i MAPPED_SAMPLE_DIR` | Path to the mapped sample directory. | required |
| `-r REF_DIR` | Path to reference directory (Path containing *_reference.json) | required |
| `-d MAX_DEPTH` | Maximum read depth taken in account at each position. | `100` |
| `-p MIN_DEPTH` | Minimum read depth of a site: sites covered by fewer reads are masked in the consensus. Values should be comprised between 1 and the maximum depth (`-d`). | `3` |
| `-s MIN_SNP_DEPTH` | Minimum read depth of a site to call a variant. Values should be comprised between 1 and the maximum depth (`-d`), and be >= the minimum depth (`-p`). | `3` |
| `-f MIN_FREQUENCY` | Minimum allele frequency to call a variant. | `0.1` |
| `-l PLOIDY` | Ploidy given to the variant caller. | `1` |
| `-m MIN_MSP_COVERAGE` | Minimum percentage (1 to 100) of the signature genes of an MSP passing `-c` to build its strain. | `80` |
| `-c MIN_GENE_COVERAGE` | Minimum fraction of a signature gene covered by reads, from 0 to 1. | `0.5` |
| `--core_size CORE_SIZE` | Number of signature genes per species (MSP) used to estimate their respective abundance. | `100` |
| `-o STRAIN_DIR` | Path to output directory. | required |
| `--kc` | Keep the consensus of the marker genes, needed to run meteor strain again on the sample. | off |
| `--tmp TMP_PATH` | Path to the directory where temporary files are stored | none |
| `-t THREADS` | Number of threads to perform variant calling. | `1` |

### meteor tree

| Option | Description | Default |
|---|---|---|
| `-i STRAIN_DIR` | Path to the strain directory. | required |
| `-g MAX_GAP` | Removes sites constitued of >= cutoff gap character. | `1.0` |
| `-s MIN_INFO_SITES` | Minimum number of informative sites in the alignment. | `10` |
| `-r` | Compute distances with the TN93 model (faster) instead of GTR. | GTR |
| `--treeshrink` | Flag the samples on abnormally long branches of each strain tree (per-gene test of TreeShrink): writes <msp>_qcflags.tsv and <msp>_shrink.tsv. | off |
| `-l SAMPLE_FILE` | Build the trees with only the samples of this file: one sample name per line (first column if tab-separated; empty lines and lines starting with # ignored). | none |
| `--groups GROUPS_FILE` | Tab-separated file of the group of each sample (sample, group), e.g. household, mother-infant pair or subject. Samples of the same group are related: strain identity thresholds are set from the unrelated pairs and the transmissibility of each MSP from the related pairs. | none |
| `--min_overlap MIN_OVERLAP` | Compare two samples on an MSP when the sites informative in both cover at least this % of its alignment. | `10.0` |
| `--min_pairs MIN_PAIRS` | Minimum number of compared sample pairs to report the sharing rate and transmissibility of an MSP. | `10` |
| `-f FORMAT` | Output image format of the trees. Choices: `png`, `svg`, `pdf`, `txt`. | no image |
| `-w WIDTH` | Output image width. | `500` |
| `-H HEIGHT` | Output image height. | `500` |
| `-o TREE_DIR` | Path to output directory. | required |
| `--tmp TMP_PATH` | Path to the directory where temporary files are stored | none |
| `-t THREADS` | Number of threads when infering each tree. | `1` |

## Citing Meteor2

Please cite the following publication if you use Meteor2:  
[Accurate profiling of microbial communities for shotgun metagenomic sequencing with Meteor2.](https://doi.org/10.1186/s40168-025-02249-w)   
Amine Ghozlane, Florence Thirion, Florian Plaza Oñate, Franck Gauthier, Emmanuelle Le Chatelier, Anita Annamalé, Mathieu Almeida, Stanislav D. Ehrlich, Nicolas Pons.
