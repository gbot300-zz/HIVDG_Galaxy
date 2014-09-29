#!/usr/bin/python

'''

Author: Gordon

11 Dec 2013

Wrapper for PEAR Illumina read merger tool.

'''

import os
import subprocess

#import utilities as util

import logging
logging.basicConfig(format= '%(asctime)-15s >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

import argparse
parser = argparse.ArgumentParser(description = ' PEAR runner logging')

## Required input
parser.add_argument('--f', help='file with forward reads')
parser.add_argument('--r', help='file with reverse reads')
parser.add_argument('--v', help='minimum-overlap')
parser.add_argument('--m', help='max-assembly-length')
parser.add_argument('--n', help='min-assembly-length')
parser.add_argument('--t', help='min-trim-length')
parser.add_argument('--q', help='quality-threshold')

args = parser.parse_args()
dargs =  vars(args)
print dargs

########################################################################################

## required inputs
cmd = 'pear -f %s -r %s -o out ' % (dargs['f'], dargs['r'])

## optional inputs
cmd += ' -v %s ' %  (dargs['v'])
cmd += ' -m %s ' % (dargs['m'])
cmd += ' -n %s ' %  (dargs['n'])
cmd += ' -t %s ' % (dargs['t'])
cmd += ' -q %s ' % (dargs['q'])

#cmd = '/opt/exp_soft/pear-0.8.1-src/src/pear -f %s -r %s -o %s' % (dargs['f']. dargs['r'], dargs['o'])

cmd += ' &> /dev/null '

subprocess.check_call(cmd, shell=True)
