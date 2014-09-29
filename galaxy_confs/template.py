#!/usr/bin/python

'''

Author: Gordon

Date: xxx 

'''

import os
import subprocess
import sys

import utilities as util

from maps import translator

import random

import logging
#logging = logging.getLogger('Logger')
#logging.setLevel('INFO')

import argparse
parser = argparse.ArgumentParser(description = 'Cycles through a set of merged Illumina reads and identifies Primer IDs. New PIDs are compared to an existing set and if mismatched by only one base, the sequence is added to the existing PID')

## TODO
	## 1. INPUTS
	## 2. FUNCTIONS
	## 3. Outputs

## Required input
parser.add_argument('--', help='')
parser.add_argument('-', action="store_true", help="")

args = parser.parse_args()
dargs =  vars(args)
#print dargs

logFilename = "log-" + str(random.randint(0,10000))
logging.basicConfig(format= '%(asctime)-15s || %(levelname)s  >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', \
filename=logFilename, filemode='w', level=logging.INFO)

def ...
	
#######################################################################################3
logging.info("")
logging.info("done...")



logging.info('moving LOG file')
subprocess.check_call('mv %s log.txt;' % (logFilename), shell=True)
logging.info('done...')










