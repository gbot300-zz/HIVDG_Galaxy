#!/usr/bin/python

'''
Wrapper for the Broad Tool AV454 which assembles reads specifically for viral genomes

Author: gordon
Date: 10 July 2013
'''

import os










'''
Script Usage

Once Arachne is fully installed and configured, you can run the module for AssembleViral454 directly. It is located in the 
src/assemble folder in the main Arachne package.

Command line:
AssembleViral454 ORG=<OrganismName> GSIZE=<GenomeSize> PRE=<PRE_DIRECTORY>  DATA=<DATA_DIRECTORY> 
			RUN=<RUN_DIRECTORY> SFF_DIR=<SFF_DIRECTORY>

To obtain a fasta of the assembly, you can use the following module (also installed with Arachne):

MakeReadableOutput PRE=<PRE_DIRECTORY> DATA=<DATA_DIRECTORY> RUN<RUN_DIRECTORY>/work/ SUBDIR=post_final

For details on what the PRE, DATA and RUN directories are, see Arachne documentation. 

The SFF directory must contain all SFF files that you want assembled, it will use all the files in the folder for the assembly.

It should create the file “assembly.bases.gz” that you can then gunzip.
'''
