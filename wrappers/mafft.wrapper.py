#!/usr/bin/python

'''
description

'''

import os
import subprocess
import argparse
import sys

## Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for mafft aligner')

parser.add_argument('--i', help='Path to FASTA input file')

args = parser.parse_args()
dargs =  vars(args)

## tool execution
cmd = ' mafft --auto --quiet %s > mafft_output ' % (dargs['i'])

print cmd
subprocess.check_call(cmd, shell=True)
