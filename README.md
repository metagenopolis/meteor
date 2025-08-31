# Meteor

[![install with conda](https://img.shields.io/conda/vn/bioconda/meteor?color=green&label=bioconda%2Fmeteor&logo=anaconda)](https://anaconda.org/bioconda/meteor)
![Anaconda downloads](https://anaconda.org/bioconda/meteor/badges/downloads.svg)
[![PyPI](https://img.shields.io/pypi/v/METEOR?label=pypi%20package)](https://pypi.org/project/meteor/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/meteor?label=meteor%20on%20PyPi)
[![dockerhub](https://img.shields.io/docker/v/aghozlane/meteor?label=aghozlane/meteor&logo=docker)](https://hub.docker.com/r/aghozlane/meteor/)
![Github Actions](https://github.com/metagenopolis/meteor/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/metagenopolis/meteor/graph/badge.svg?token=AXAEIUY7DX)](https://codecov.io/gh/metagenopolis/meteor)
[![DOI](https://zenodo.org/badge/722959292.svg)](https://zenodo.org/doi/10.5281/zenodo.10912587)

## Introduction

Meteor is a plateform for quantitative metagenomics profiling of complex ecosystems.
Meteor relies on genes catalogue to perform species-level taxonomic profiling (Bacteria, Archaea and Eukaryotes), functional analysis and strain-level population structure inference.


## Dependencies

Besides python packages dependencies, Meteor requires:
- python 3.10+
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

## Citing Meteor2

Please cite the following publication if you use Meteor2:  
[Accurate profiling of microbial communities for shotgun metagenomic sequencing with Meteor2.](https://doi.org/10.21203/rs.3.rs-6122276/v1)   
Amine Ghozlane, Florence Thirion, Florian Plaza Oñate, Franck Gauthier, Emmanuelle Le Chatelier, Anita Annamalé, Mathieu Almeida, Stanislav D. Ehrlich, Nicolas Pons.
