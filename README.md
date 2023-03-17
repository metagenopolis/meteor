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
- g++ and make
- [bowtie2](https://github.com/BenLangmead/bowtie2)

## Installation

Meteor is available with bioconda which includes all its dependencies:
```
conda install -c bioconda meteor
```

And with pip:
```
pip install meteor
```
You can test the installation of meteor with:
```
meteor test
```
## Getting started

### Download a reference catalogue

Meteor requires to download locally a microbial gene catalogue. Several catalogues are currently available:

|  Microbial gene catalogue | \<name\> | Genes count | Metagenomic Species Pan-genomes (MSPs) |Size (Gb)  | Description  |
|:---|:---:|:---:|:---:|:---:|:---:|
|  *Gallus gallus domesticus* | chicken_caecal  | 13.6M  | 2420 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FHPJH5)
| *Homo sapiens oral*  |  human_oral | 8.4M  | 853 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/WQ4UTV)
| *Homo sapiens gut* |  human_gut | 10.4M  | 1990 | |[link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FLANUP)
|  *Rattus norvegicus* | rat_gut  | 5.9M | 1627 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/GVL2EE)
|  *Mus musculus*  | mouse_gut  | 5.0M  | 1252 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/L11MXM)
|  *Sus domesticus* | pig_gut  | 9.3M  | 1523 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/OPAULL)
|  *Oryctolagus cuniculus* | rabbit_gut  | 5.7M | 1053 | | [link](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/5EJKAS)

These references can be downloaded with following command:
```
meteor download -i <name> -c -o <refdir>
```

### Import custom reference

Users can also import custom gene catalogue with the command:
```
meteor build -i <fastafile> -n <name> -o <refdir> -t <threads>
```

### Import fastq

Meteor requires a first of fastq indexing:
```
meteor fastq -i <fastqdir>  [-p paired reads] -n <projectname> -o <outputdir>
```
When multiple sequencing were performed for a library, the option -m allows to group these samples.
Example:

Illumina_lib1-**SAMPLE_01**.fastq
Illumina_lib1-**SAMPLE_02**.fastq
Illumina_lib2-**SAMPLE_01**.fastq
Illumina_lib2-**SAMPLE_02**.fastq

In this case, the command would be:
```
meteor fastq -i ./  -m SAMPLE_\\d+ -n projectname -o outputdir
```
Which group samples the same id.

### Profiling

The profiling can be obtained from read:
```
meteor mapping -i <fastqdir> -r <refdir> -o <mappingdir>
```

### Taxonomic analysis
The taxonomic profiling can be performed several level of accuracy (gene, MSP, SuperKingdom, Phylum,Class, Order, Family, Genus, Specie)
```
meteor profiler -i <mappingdir> -l <accuracy> -o <countingdir>
```

### Functionnal analysis

Meteor provides an ARDs annotation based on (MUSTARD)[https://www.nature.com/articles/s41564-018-0292-6] and (KEGG)[https://academic.oup.com/nar/article/36/suppl_1/D480/2507484]
```
meteor profiler -i <mappingdir> -o <countingdir>
```

## The METEOR team
The main contributors to METEOR:
* Nicolas Pons
* Franck Gauthier
* Amine Ghozlane
* Florence Thirion
* Florian Plaza Oñate

## Acknowledgements
Special thanks to the following people:
* Emmanuelle Le Chatelier
* Mathieu Almeida