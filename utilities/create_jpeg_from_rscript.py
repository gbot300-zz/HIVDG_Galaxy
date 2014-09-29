#!/usr/bin/python

''' 
Author: Gordon

Date: 07 Oct 2013

Script to to R-file as input and create a jpeg instead of a PDF image. Originally created to redo the Heatmaps from the 454 pipeline

'''

import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description = 'Create Jpeg versions of PDF image files')

## Required input
parser.add_argument('--rfile')
parser.add_argument('--dfile')

args = parser.parse_args()
dargs =  vars(args)
print dargs

with open(dargs['rfile'], 'r') as f:
	lines = f.readlines()

with open('Rscript_jpeg.r', 'w') as w:
	for line in lines:
		if line.startswith('pdf('):
			w.writelines(['jpeg(filename = "jpg_plot.jpg", quality=100)', os.linesep])
		elif line.startswith('pmat <-read.delim('):
			x = 'pmat <-read.delim("%s",header=TRUE,row.names=1, skip=1)' % (dargs['dfile'])
			w.writelines([x, os.linesep])
		else:
			w.write(line)


cmd = ' R CMD BATCH Rscript_jpeg.r'

cmd += ' &> /dev/null '

print cmd

subprocess.check_call(cmd, shell=True)
