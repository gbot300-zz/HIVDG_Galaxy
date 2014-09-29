#!/usr/bin/python

'''
Script to split FastQ files into a FASTA reads file and a QUAL quality file

Author: Gordon
Date: 23 July 2013

Dependencies: BioPython

'''
import os
import argparse
from Bio import SeqIO


def reWrite(fx, wx, space=False):
	if space:
		s=' '
	else:
		s=''
	with open(fx, 'r') as f:
		lines = f.readlines()
		lines2 = []
		c = 0
		for line in lines:
			c +=1
			if not line.startswith('>'):
				lines2.append(line.strip(os.linesep) +s )
			# had to add this stupid hack because Mosaik doesnt pick up a file sa=tarting with an empty line as a fasta file...
			elif c != 1:
				zz = os.linesep+'>'
				lines2.append(line.replace('>', zz))
			else:
				lines2.append(line)

	with open(wx, 'w') as w:
		w.writelines(lines2)

# Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Splits FastQ file into FASTA reads file and QUAL quality file')
parser.add_argument('-i',      	help = 'FastQ input file')
parser.add_argument('of', 	help = 'Fasta out file')
parser.add_argument('qf',       help = 'Qual out file')
args = parser.parse_args()
dargs =  vars(args)

#print dargs

basename 	= os.path.splitext(dargs['i'])[0]
fqpath		= os.path.join(os.getcwd(), dargs['i'])
#faname		= basename + ".fasta"
#quname		= basename + ".qual"	
#fapath		= os.path.join(os.getcwd(), faname)
#qupath		= os.path.join(os.getcwd(), quname)

SeqIO.convert(fqpath, "fastq", 'xxx.fasta.tmp', "fasta")
SeqIO.convert(fqpath, "fastq", 'yyy.qual.tmp', "qual")

# Added this to cut the newlines characters out of sequences which, for some dumb reason, SeqIO adds 
reWrite('xxx.fasta.tmp', 'xxx.fasta')
reWrite('yyy.qual.tmp', 'yyy.qual', space=True)

#os.rename(fapath, dargs['of'])
#os.rename(qupath, dargs['qf'])
