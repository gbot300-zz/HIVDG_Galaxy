#!/usr/bin/python

'''

Author: Gordon

Date 11 Dec 2013

'''

import os
import subprocess

import dendropy
from dendropy import treecalc

import argparse
parser = argparse.ArgumentParser(description = 'Calculates euclidian between reference sequence and test sequence in a tree')

## Required input
parser.add_argument('--f', help='Alignment file, with reference sequence as first sequence')
parser.add_argument('--r', help='Reference sequence name')
parser.add_argument('--p', help='protein of nucleotide')

args = parser.parse_args()
dargs =  vars(args)
print dargs

#####################################################################################

ref = dargs['r']

if dargs['p'] == "protein":
	cmd = ' FastTreeMP < %s > out.tree -quiet' % (dargs['f'])
else:
	cmd = ' FastTreeMP -nt < %s > out.tree -quiet' % (dargs['f'])

subprocess.check_call(cmd, shell=True)

tree = dendropy.Tree.get_from_path('out.tree', 'newick')

pdm = treecalc.PatristicDistanceMatrix(tree)

with open('out.csv', 'w') as w:
	for node in tree.taxon_set:
		if node.label.replace(' ', '_') == ref:
			w.write('seq, score with Ref: %s\n' % ref)
			for node2 in tree.taxon_set:
				if not node2.label.replace(' ', '_') == ref:
					w.write(' %s, %2f \n' % (node2.label.replace(' ', '_'), pdm(node, node2))) # fix this part 
	

















