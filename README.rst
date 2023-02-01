=======
Meteor
=======
## Notes for Meteor Version 3.2

Introduction
============

Meteor is a plateform for quantitative metagenomics profiling of complex ecosystems.
Meteor relies on genes catalogue to perform specie level taxonomic assignments and functional analysis. 

Installation
============

You can install Meteor using the pip command::

    pip install meteor

Meteor is also available with bioconda::

    conda install -c bioconda meteor

Getting started
===============

Download reference
------------------

Several gut microbial gene catalogue were designed.

+-------+------------+
| Name  | Reference  |
+=======+============+
| human |  ref  |
+-------+------------+
| mouse |  ref  |
+-------+------------+
| pig |  ref  |
+-------+------------+
| chicken |  ref  |
+-------+------------+

These references can be downloaded and indexed::

    meteor download <name> <outputdir>

Import custom reference
-----------------------

Users can also import custom gene catalogue with the command::

    meteor reference ....

Import fastq
------------

Meteor requires a first of fastq indexing::

    meteor fastq -i <name>  -p ....


Profiling
----------

The taxonomic profiling can be obtained from read::

    meteor mapping ..


Analysis
========