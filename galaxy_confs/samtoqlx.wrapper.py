#!/usr/bin/python

'''
description

'''

import os
import subprocess
import argparse

## Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for samToQlx.pl')

## required args
parser.add_argument('--sam', help='SAM file to be converted')
parser.add_argument('--reads', help="FASTA format reads file")
parser.add_argument('--quals', help="QUAL format quality file")
parser.add_argument('--ref', help="FASTA format reference file")

args = parser.parse_args()
dargs =  vars(args)
print dargs

## tool execution

cmd = '/usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/samToQlx.pl '

## required inputs
cmd += ' %s %s %s %s out.qlx ' % (dargs['sam'], dargs['reads'], dargs['quals'], dargs['ref'])

print cmd
subprocess.check_call(cmd, shell=True)

