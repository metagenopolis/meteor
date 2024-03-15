# Meteor

[![install with conda](https://img.shields.io/conda/vn/aghozlane/meteor?color=green&label=aghozlane%2Fmeteor&logo=anaconda)](https://anaconda.org/aghozlane/meteor)
[![PyPI](https://img.shields.io/pypi/v/METEOR?label=pypi%20package)](https://pypi.org/project/meteor/)
[![dockerhub](https://img.shields.io/docker/v/aghozlane/meteor?label=aghozlane/meteor&logo=docker)](https://hub.docker.com/r/aghozlane/meteor/)
![Github Actions](https://github.com/metagenopolis/meteor/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/metagenopolis/meteor/graph/badge.svg?token=AXAEIUY7DX)](https://codecov.io/gh/metagenopolis/meteor)

## Introduction

Meteor is a plateform for quantitative metagenomics profiling of complex ecosystems.
Meteor relies on genes catalogue to perform species-level taxonomic profiling, functional analysis and strain-level population structure inference.


## Dependencies

Besides python packages dependencies, Meteor requires:
- python 3.10+
- [bowtie2 2.3.5+](https://github.com/BenLangmead/bowtie2)
- [bcftools 0.1.19+](https://samtools.github.io/bcftools/)
- [bedtools 2.18.0+](https://bedtools.readthedocs.io/en/latest/index.html)
- [FastTree 1.9.0+](http://www.microbesonline.org/fasttree/)

## Installation

Meteor is available with conda which includes all its dependencies:
```
conda create --name meteor -c conda-forge -c bioconda  -c aghozlane meteor
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

### 1. Download or build a custom reference
-------------------------------------------

Meteor requires to download locally a microbial gene catalogue specif, either in 'full' or 'light' version. The 'full' version contains all genes of the catalogue, whereas the 'light' version contains only the marker genes that will be used to infer species abundance profiles. Of note, no functional profiling can be performed when using the 'light' version of a catalogue.

Ten catalogues are currently available:

|  Microbial gene catalogue | \<name\> | Genes count (M) | Metagenomic Species Pan-genomes (MSPs) |Size (full) (GB) | Size (light) (GB)  | Description  |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
|  *Felis catus* | cat_gut  | 1.3  | 344 | 2.7 | 0.9 |[link](https://zenodo.org/records/10719585)
|  *Gallus gallus domesticus* | chicken_caecal  | 13.6  | 2420 | 22.8 | 4.8 |[link](https://zenodo.org/records/10719564)
|  *Canis lupus familiaris* | dog_gut  | 0.95  | 234 | 1.7 | 0.3 |[link](https://zenodo.org/records/10719585)
| *Homo sapiens gut* |  human_gut | 10.4  | 1990 | 15.1 | 3.2 |[link](https://zenodo.org/records/10719553)
| *Homo sapiens oral*  |  human_oral | 8.4  | 853 | 16.4 | 3.2 |[link](https://zenodo.org/records/10719586)
| *Homo sapiens skin*  |  human_skin | 2.9  | 392 | 4.6 | 0.9 |[link](https://zenodo.org/records/10719613)
| *Mus musculus*  | mouse_gut  | 5.0  | 1252 | 10.3 | 3.4 |[link](https://zenodo.org/records/10719617)
| *Oryctolagus cuniculus* | rabbit_gut  | 5.7 | 1053 | 11.1 | 3.5 |[link](https://zenodo.org/records/10719606)
| *Rattus norvegicus* | rat_gut  | 5.9 | 1627 | 8.4 | 2.0 |[link](https://zenodo.org/records/10719596)
| *Sus domesticus* | pig_gut  | 9.3  | 1523 | 8.4 | 378 |[link](https://zenodo.org/records/10719591)

These references can be downloaded with the following command:
```
meteor download -i <name> -c -o <refdir>
```
The 'light' catalogues are available with the tag (--fast) :
```
meteor download -i <name> -c --fast -o <refdir>
```

Users can also import custom gene catalogue with the command:
```
meteor build -i <fastafile> -n <name> -o <refdir> -t <threads>
```

### 2. Import fastq
-------------------
Meteor requires a first of fastq indexing:
```
meteor fastq -i <fastqdir>  [-p paired reads] -n <projectname> -o <outputdir>
```
When multiple sequencing are available for a library, the option -m allows to group these samples.
Example:

Illumina_lib1-**SAMPLE_01**.fastq <br />
Illumina_lib1-**SAMPLE_02**.fastq <br />
Illumina_lib2-**SAMPLE_01**.fastq <br />
Illumina_lib2-**SAMPLE_02**.fastq <br />

In this case, the following command will group these samples the same library:
```
meteor fastq -i ./  -m SAMPLE_\\d+ -n projectname -o outputdir
```

### 3. Mapping
----------------
The raw fastq files are mapped against a catalogue to generate a gene count table with the following command:
```
meteor mapping -i <fastqdir/sampledir> -r <refdir> -o <mappingdir>
```
We recommend to first filter out reads with low-quality, length < 60nt or belonging to the host.

### 4. Taxonomic and functional profiling
-------------------------

Genes from the catalogue are clustered into Metagenomic Species Pangeomes (MSP) with [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712), and are functionnaly annotated against [KEGG r107](https://academic.oup.com/nar/article/36/suppl_1/D480/2507484), [DBcan](https://academic.oup.com/nar/article/51/W1/W115/7147496?login=true) (carbohydate active enzymes) and [MUSTARD](https://www.nature.com/articles/s41564-018-0292-6) (antibiotic resistant determinants).

 MSP and functional profiles are computed from the gene count table with the following command:

```
meteor profile -i <mappingdir> -o <profiledir> -r <refdir> -n coverage
```

The "-n" parameter ensures read count normalization for gene length. If omitted, no normalization will be performed on the gene table.

This profiling step will generate:
- species abundance table;
- ARD abundance table (full catalogue only);
- DBCAN abundance table (full catalogue only);
- Gut Metabolic Modules ([GMM](https://www.nature.com/articles/nmicrobiol201688)) abundance table (from the KO annotation) (full catalogue only).


### 5. Merging

To merge output from different samples into a single table, use the following command:

```
meteor merge -i <profiledir> -o <mergingdir> --fast
```

The '--fast' parameter prevent merging of the gene count tables, so that only species and functions table will be merged.

### 5. Strain profiling
-------------------------


## The METEOR team
The main contributors to METEOR:

* Franck Gauthier
* Amine Ghozlane
* Florian Plaza Oñate
* Nicolas Pons
* Florence Thirion


## Acknowledgements
Special thanks to the following people:
* Mathieu Almeida
* Emmanuelle Le Chatelier
