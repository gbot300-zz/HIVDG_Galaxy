#!/usr/bin/python

'''
description

'''

import os
import subprocess
import argparse

## Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for bin_by_seqstring.py')

## required args
parser.add_argument('--fa', help='Merged reads in fasta format')

## optional args
parser.add_argument('--consensus', action="store_true", help='Do you want to find a consensus sequence for the sequences in each bin?')

args = parser.parse_args()
dargs =  vars(args)
print dargs
## tool execution
cmd = '/usr/bin/python /opt/exp_soft/galaxy-tools/custom_tools/XXX '

## required inputs
if dargs['fa']:
	cmd += ' --fa %s ' % (dargs['fa'])
elif dargs['fq']:
	cmd += ' --fq %s ' % (dargs['fq'])

print cmd
subprocess.check_call(cmd, shell=True)

