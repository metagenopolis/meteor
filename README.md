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

## Getting started


### Download reference

Several gut microbial gene catalogue were designed.

These references can be downloaded and indexed:
```
meteor download <name> <outputdir>
```

### Import custom reference

Users can also import custom gene catalogue with the command:
```
meteor reference ....
```

### Import fastq

Meteor requires a first of fastq indexing:
```
meteor fastq -i <name>  -p ....
```

### Profiling

The taxonomic profiling can be obtained from read::
```
meteor mapping ..
```

### Analysis

