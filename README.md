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
- python>3.11.0
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

## Installation

Meteor is available with bioconda which includes all its dependencies:
```
conda install -c bioconda meteor
```

And with pip with a recent Python >= 3.11:
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
2. **Import the raw fastq files**
3. **Map reads against the reference catalogue**
4. **Profile taxonomical or functionnal abundances**


### 1. Download or build a custom reference
-------------------------------------------

Meteor requires to download locally a microbial gene catalogue. Several catalogues are currently available:

|  Microbial gene catalogue | \<name\> | Genes count (M) | Metagenomic Species Pan-genomes (MSPs) |Size (GB) | Taxonomy catalogue size (MB)  | Description  |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
|  *Gallus gallus domesticus* | chicken_caecal  | 13.6  | 2420 | 15.5 | 628 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FHPJH5)
| *Homo sapiens oral*  |  human_oral | 8.4  | 853 | 10.1 | 179 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/WQ4UTV)
| *Homo sapiens gut* |  human_gut | 10.4  | 1990 | 9.3 | 391 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FLANUP)
|  *Mus musculus*  | mouse_gut  | 5.0  | 1252 | 5.5 | 347 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/L11MXM)
|  *Oryctolagus cuniculus* | rabbit_gut  | 5.7 | 1053 | 5.9 | 199 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/5EJKAS)
|  *Rattus norvegicus* | rat_gut  | 5.9 | 1627 | 5.1 | 348 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/GVL2EE)
|  *Sus domesticus* | pig_gut  | 9.3  | 1523 | 8.4 | 378 |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/OPAULL)


These references can be downloaded with following command:
```
meteor download -i <name> -c -o <refdir>
```
We also created smaller catalogue designed exclusively for taxonomical profiling. They are available with the tag (-t) :
```
meteor download -i <name> -c -t -o <refdir>
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

### 4. Abundance profiling
-------------------------

####  **Taxonomical analysis**

The taxonomical profiling can be performed at several level of accuracy (gene, MSP, SuperKingdom, Phylum, Class, Order, Family, Genus, Specie). MSPs were calculated with [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712).
Their abundances can be profiled with the following command:
```
meteor profiler -i <mappingdir> -l <accuracy> -o <countingdir>
```

#### **Functionnal analysis**

Meteor provides an ARDs annotation based on [MUSTARD](https://www.nature.com/articles/s41564-018-0292-6) and a pathway annotation based on [KEGG r107](https://academic.oup.com/nar/article/36/suppl_1/D480/2507484). Their abundances can be profiled with the following command:
```
meteor profiler -i <mappingdir> -a <annotation> -o <countingdir>
```

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
