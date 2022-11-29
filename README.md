=======================================================================
Release notes for Meteor Version 3.2
=======================================================================

Files and directories
=====================
In the package, files are organized as follows:

  -+- README:    This file
   |
   +- Quick_start_Meteor.pdf :   Meteor Usage instructions
   |
   +- COPYING:   License
   |
   +- AUTHORS:   List of authors
   |
   +- data_preparation_tools/ --> programs dedicated to data preparation  
   |
   +- meteor-pipeline/ --> sources and other files to build the Meteor pipeline.
      |
      +- src/ --> c++ sources and makefile.
      |
      +- meteor.rb:   ruby source file.


======================================================================
Installation notes for Meteor Version 3.2
======================================================================

 This file describes the building and installation instructions from the SOURCES for UNIX like systems.
 NB: Building and running Meteor on MS Windows system has not been tested,
     but should easily be achieved as Meteor uses standard C++ code without any exotic dependencies.

Supported Systems
=================

 Will work in most modern GNU/Linux, and macOS (10.5 or later with Intel processor) systems.

Building and running pre-requisite
==================================

 1. Compiler collections with c++98 support (e.g. g++ or clang++)

 2. Ruby interpreter version 1.9 or later, with "inifile" rubygem module version 3.0.0 or later

 3. bowtie2 version 2.3.5.1 or later, and bowtie 1.2.3 or later if you plan to use color space SOLiD data.

Installation
============

 On a GNU/Linux system, we recommand to have g++, make and ruby installed via the package manager of your
 distribution.
 e.g.:

    on Debian/Ubuntu plateforms (as root or sudoer):
    apt-get install make g++ ruby

    on Redhat/Fedora/CentOS plateforms (as root or sudoer):
    yum install make gcc-c++ ruby

 - To build Meteor, open a terminal and execute the following commands:

    cd meteor-pipeline/src/
    make

 This will create two binaries (meteor-counter and meteor-profiler) in meteor-pipeline/src/build/

 - Now copy the meteor pipeline executables in a folder present in the PATH environment variable.

    e.g. as root or sudoer (this will install meteor for all user):
    cp build/meteor-counter build/meteor-profiler ../meteor.rb /usr/local/bin

    Or as common user (every user of meteor will have to do this):
    [ ! -d ~/bin ] && mkdir ~/bin
    cp build/meteor-counter build/meteor-profiler ../meteor.rb ~/bin

 - Then install the ruby module inifile (>= 3.0.0).

    As root or sudoer:
    gem install inifile

    As common user (every user of meteor will have to do this):
    gem install inifile --user-install


======================================================================
Usage of Meteor Version 3.2
======================================================================

Data preparation
================

Use the program MeteorReferenceBuilder.rb to create the bowtie index of the reference gene catalog(s)
Then organize the sample data with MeteorImportFastq.rb or MeteorImportCSFastaQual.rb

These programs are located in the folder data_preparation_tools/

Using Meteor 
============

   meteor.rb [list_of_arguments]
   meteor-profiler [list_of_arguments]

meteor.rb drives the mapping and counting processes for one sample.
So launch as many instance of meteor.rb as the number of samples in your project.
Then use meteor-profiler to gather the counting results of all samples in a unique matrix.

NOTE: It is necessary to include the executable directory in the PATH environment variable.

See Meteor quick start documentation for more details.

