#!/usr/bin/python

'''

Author: Gordon GBOT300@GMAIL.COM	

Date: 9 July 2014

Script reads a file with sequences into an ordered dictionary and removes all sequences that contain degenerate bases

'''

import os
import logging
import argparse
import random

parser = argparse.ArgumentParser(description = 'Reads through a file and removes all sequences that contain degenerate bases. Outputs remaining sequences to file')

parser.add_argument('--fa', help='Reads in fasta format')
parser.add_argument('--fq', help='Reads in fastq format')

args = parser.parse_args()
dargs =  vars(args)
#print dargs

logFilename = "log_degenremover-" + str(random.randint(0,10000))
logging.basicConfig(format= '%(asctime)-15s || %(levelname)s  >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', \
		filename=logFilename, filemode='w', level=logging.INFO)

import utilities as util

if dargs['fa']:
	fname		= dargs['fa']
	order, seqs	= util.clean_seqs(fname) 
else:
	fname				= dargs['fq']
	order, seqs, quals  = util.clean_seqs(fname, fastq=True)

new_order = order[:]

for k in order:
	seq = seqs[k]
	for db in ['M', 'R', 'W', 'S', 'Y', 'K', 'B', 'D', 'H', 'V', 'N']:
		if (db in seq) and (k in new_order):
			new_order.remove(k)

if dargs['fa']:
	util.create_fasta(new_order, seqs, 'fasta_degen_removed')
else:
	util.create_fastq(new_order, seqs, quals, 'fastq_degen_removed')

