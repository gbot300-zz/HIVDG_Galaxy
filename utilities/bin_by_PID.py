#!/usr/bin/python

'''

Author: Gordon

Date 13 Dec 2013

'''

import os
import subprocess

import utilities as util

from fuzzywuzzy import fuzz, process

import argparse
parser = argparse.ArgumentParser(description = 'Cycles through a set of merged Illumina reads and identifies Primer IDs. New PIDs are compared to an existing set and if mismatched by only one base, the sequence is added to the existing PID')

## TODO if the first primer found has a mismatch, all subsequent "mismatches" to that primer will have their sequences re-assigned to the 
## mismatched primer. Think of a way to fix this
## TODO rewrite your scripts so that importing them into another one is possible. such as for the fastq_to_fasta_qual.py one
## TODO test similar PIDs by pairwise alignments with exisitng files

## Required input
parser.add_argument('--f', help='Merged reads in fastq format')
#parser.add_argument('--r', help='Reference sequence name')

args = parser.parse_args()
dargs =  vars(args)
#print dargs

#######################################################################################3

rev_uni_primer 	= 'CTGAGCGTGTGGCAAGGC'

snames, seqs, quals 	= util.clean_seqs(os.path.join(os.getcwd(), dargs['f']), fastq=True)
#print snames

unique_primers 	= []
no_primers	= []

bins 		= {}

for sname in seqs.keys():
	#print sname
	if fuzz.partial_ratio(rev_uni_primer, seqs[sname]) > 80:
		primer = seqs[sname][-27:-18]
		x = process.extract(primer, unique_primers, limit=1)
		if not x == [] and x[0][1] > 87:
			bins[x[0][0]].append(sname)
		else:
			unique_primers.append(primer)
			bins[primer] = [sname]

	else:
		no_primers.append(sname)

#print 'nop: ', no_primers
#print 'bins: ', bins

for key in bins.keys():
	odir = os.path.join(os.getcwd(), 'bins', key[0], key[1], key[2], key[3], key[4], key[0], key[5], key[6], key[7], key[8])
	util.ensure_dir_exists(odir)
	wpath = os.path.join(odir, key + '_sequences.fastq')
	util.create_fastq(bins[key], seqs, quals, wpath)

for seq in no_primers:
	odir = os.path.join(os.getcwd(), 'bins', 'no_primers')
	util.ensure_dir_exists(odir)
	wpath = os.path.join(odir, 'no_primers_sequences.fastq')
	util.create_fastq(no_primers, seqs, quals, wpath)

cmd = ' zip -r bins.zip bins  &> /dev/null;  '

#'vprofiler.pl -i input.txt -o Vpro -noendvariant=10 -nt -codon -haplo -haploseq '
subprocess.check_call(cmd, shell=True)









