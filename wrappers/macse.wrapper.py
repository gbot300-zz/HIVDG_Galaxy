#!/usr/bin/python

'''

Author: Gordon

Date: Aug 2014

Wrapper for MACSE codon alignment tool; http://bioweb.supagro.inra.fr/macse/index.php?menu=intro&option=intro

'''


import os
import subprocess
import sys
import random

import logging
#logging = logging.getLogger('Logger')
#logging.setLevel('INFO')

import argparse
parser = argparse.ArgumentParser(description = 'Wrapper for MACSE Codon Alignment tool')

## TODO

## Required input
parser.add_argument('--seq', help='input sequence in FASTA format')

## Optional Input
parser.add_argument('--ext_gap_ratio', help='the gap cost is multiplied by this ratio (< 1) to define the cost of gaps at the begining/end of sequences.')
parser.add_argument('--fs', help=' cost of a frameshift within the sequ`ence')
parser.add_argument('--gap_ext', help=' cost of a gap extension')
parser.add_argument('--gap_opt', help=' cost of creating a gap')
parser.add_argument('--stop', help=' cost of a stop codon not at the end of the sequence')

args = parser.parse_args()
dargs =  vars(args)
#print dargs

logFilename = "log-macse" + str(random.randint(0,10000))
logging.basicConfig(format= '%(asctime)-15s || %(levelname)s  >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', \
filename=logFilename, filemode='w', level=logging.INFO)
#######################################################################################

cmd = 'java -jar /opt/exp_soft/galaxy-tools/custom_tools/bin/macse_v1.01b.jar -prog alignSequences -out_AA macse_out_AA.fasta -out_NT macse_out_NT.fasta'
for a in dargs:
	if dargs[a]:
		cmd += ' -%s %s ' % (a, dargs[a]) 

#cmd += '  > /dev/null '
print cmd
subprocess.check_call(cmd, shell=True)

		

#######################################################################################
logging.info("")
logging.info("done...")



logging.info('moving LOG file')
subprocess.check_call('mv %s log.txt;' % (logFilename), shell=True)
logging.info('done...')










