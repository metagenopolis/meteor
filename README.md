# Meteor

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/metaphlan/README.html)
[![Latest Release](https://forgemia.inra.fr/metagenopolis/meteor/-/badges/release.svg)](https://forgemia.inra.fr/metagenopolis/meteor/-/releases)
[![pipeline status](https://forgemia.inra.fr/metagenopolis/meteor/badges/dev/pipeline.svg)](https://forgemia.inra.fr/metagenopolis/meteor/-/commits/dev)
[![coverage report](https://forgemia.inra.fr/metagenopolis/meteor/badges/dev/coverage.svg)](https://forgemia.inra.fr/metagenopolis/meteor/-/commits/dev)
[![pylint](https://forgemia.inra.fr/metagenopolis/meteor/-/jobs/artifacts/dev/raw/pylint/pylint.svg?job=pylint)](https://forgemia.inra.fr/metagenopolis/meteor/-/jobs/artifacts/dev/raw/pylint/pylint.log?job=pylint)

## Introduction

Meteor is a plateform for quantitative metagenomics profiling of complex ecosystems.
Meteor relies on genes catalogue to perform specie level taxonomic assignments and functional analysis.

Check the [wiki](https://forgemia.inra.fr/metagenopolis/meteor/-/wikis/home) for more information.
If you use meteor , please cite:



## Dependencies

Meteor requires:
- python>=3.10
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [bcftools](https://samtools.github.io/bcftools/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)
- [ete3](http://etetoolkit.org/)

## Installation

Meteor is available with conda which includes all its dependencies:
```
conda create --name meteor -c conda-forge -c bioconda  -c aghozlane meteor
```

Or with pip with a recent Python>=3.10 from the cloned directory:
```
pip install .
```
You can test the installation of meteor with:
```
meteor test
```
## Getting started

A basic usage of meteor will require to:
1. **Download or build a reference catalogue**
2. **Import the raw fastq files**
3. **Map reads against the reference catalogue**
4. **Profile taxonomical or functional abundances**
5. **Strain profiling**

### 1. Download or build a custom reference
-------------------------------------------

Meteor requires to download locally a microbial gene catalogue. Several catalogues are currently available:

|  Microbial gene catalogue | \<name\> | Genes count (M) | Metagenomic Species Pan-genomes (MSPs) |Size (GB) | Taxonomy catalogue size (GB)  | Description  |
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

These references can be downloaded with following command:
```
meteor download -i <name> -c -o <refdir>
```
We also created smaller catalogue designed exclusively for taxonomical profiling. They are available with the tag (-t) :
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
The mapping is performed with the raw fastq file against a catalogue reference with the following command:
```
meteor mapping -i <fastqdir> -r <refdir> -o <mappingdir>
```
We recommand to filter out reads of length < 60nt.

### 4. Abundance profiling
-------------------------

####  **Taxonomical analysis**

The taxonomical profiling can be performed at several level of accuracy (gene, MSP, SuperKingdom, Phylum, Class, Order, Family, Genus, Specie). MSPs were calculated with [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712).
Their abundances can be profiled with the following command:
```
meteor profile -i <mappingdir> -r <refdir> -l <accuracy> -o <profiledir>  -n coverage
```

#### **Functional analysis**

Meteor provides  a pathway annotation based on [KEGG r107](https://academic.oup.com/nar/article/36/suppl_1/D480/2507484), a CAZyme annotation based on [DBcan](https://academic.oup.com/nar/article/51/W1/W115/7147496?login=true) and an ARDs annotation based on [MUSTARD](https://www.nature.com/articles/s41564-018-0292-6). Their abundances can be profiled with the following command:
```
meteor profile -i <mappingdir> -a <annotation> -o <countingdir>
```

### 5. Strain profiling
-------------------------

#### **Functional analysis**

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
