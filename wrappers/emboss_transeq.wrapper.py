#!/usr/bin/python

'''

Author: Gordon

Date: July 2014

Wrapper for the EMBOSS tool transeq (at the time written for Emboss 6.6).

Translates DNA / RNA to AA

'''

def translate(fname):

	cmd = "transeq %s transeq_out.pep 2>&1" % (fname)
	subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":

	import subprocess

	import argparse
	parser = argparse.ArgumentParser(description = 'Translates DNA / RNA sequences to AA')

	parser.add_argument('--fa', help='Reads in fasta format')

	args = parser.parse_args()
	dargs =  vars(args)

	translate(dargs['fa'])
