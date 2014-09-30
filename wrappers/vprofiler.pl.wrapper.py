#!/usr/bin/python

'''

Python wrapper for v-profiler.pl file. The wrapper was required since the the Perl script only takes a "basename" output option and builds output files from there. I had difficulties getting Galaxy to locate the correct files

'''

import shutil
import subprocess
import argparse
import os
import ast
import cfg as cfg
#import codecs

# Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for vphaser.pl, in preparation for nucleotide, \
                                                codon and haplotype identification with V-Profiler')

## Required input
parser.add_argument('--sample_name')
parser.add_argument('--iqlx',            help = 'QLX format input file from RC454 aligner')
parser.add_argument('--icalls',          help = 'TXT format input file from V-Phaser')
parser.add_argument('--iconsen',         help = 'FASTA format input file of consensus sequence')

## Semi-required input (at least one of the options 'nt', 'codon' or 'haplo' below needs to be specified)
parser.add_argument('--nt')
parser.add_argument('--codon')
parser.add_argument('--igenes',          help = 'TXT format input file of gene positions')
parser.add_argument('--haplo')
parser.add_argument('--iregion',         help = 'TXT format input file of haplotype regions')

## optional input
#parser.add_argument('--haplosep')
parser.add_argument('--haploseq')
parser.add_argument('--cq')
parser.add_argument('--consseq')
parser.add_argument('--conshaplo')
parser.add_argument('--noendvariant')
parser.add_argument('--long_args')
#parser.add_argument('--series') #, nargs = '+')
#parser.add_argument('--sample_order') #, nargs = '+')

args = parser.parse_args()
dargs =  vars(args)
print dargs
## create input text file with a list of the input file names to adhere with structure of vprofiler.pl

#try:
#	## stick together the parts of the dictionary that have been stuck into a list with the "nargs" option in the argparse section above
#	for item in dargs['series']:
#    		x += item
#
#	## ugliest hack on earth to reform the string 'x' into a format that can be coerced into a dictionary
#	z=x.replace("{", "{'").replace("}", "'}").replace(":","':").replace(",","','").replace(":/", ":'/").replace("}',","},").replace("}'}","}}").replace("''","'")
#
#	## coerce into dictionary
#
#	series = ast.literal_eval(z)
#
#	## same with the sample_order
#	p=''.join(dargs['sample_order'])
#
#	pp = p.replace("[", "['").replace(",", "','").replace("]", "']")
#
#	sample_order = ast.literal_eval(pp)


#except TypeError:
#	series = None

with open(dargs['long_args'], 'r') as r:
	lines =  r.readlines()
	#print lines
	for i in range(len(lines)-1):
		if lines[i].strip() == '-series':
			series = ast.literal_eval(lines[i+1].strip())
		if lines[i].strip() == '-sample_order':
			sample_order = ast.literal_eval(lines[i+1].strip())
	print series
	#print sample_order

if series:
	with open('input_run.txt', 'w') as w:
		w.writelines(['>Qlx', os.linesep, 
			dargs['iqlx'], '\t', dargs['sample_name'], os.linesep])
		for item in sample_order:
			w.writelines([series[item]['iqlx'], '\t', item, os.linesep])
		w.writelines(['>VPhaser', os.linesep,
			dargs['icalls'], '\t', dargs['sample_name'], os.linesep])
		for item in sample_order:
                        w.writelines([series[item]['icalls'], '\t', item, os.linesep])
		w.writelines(['>Consensus', os.linesep,
			dargs['iconsen'], '\t', dargs['sample_name'], os.linesep])
		for item in sample_order:
                        w.writelines([series[item]['iconsen'], '\t', item, os.linesep])
else:
	with open('input_run.txt', 'w') as w:
		w.writelines(['>Qlx', os.linesep, 
				dargs['iqlx'], '\t', dargs['sample_name'], os.linesep,
				'>VPhaser', os.linesep,
				dargs['icalls'], '\t', dargs['sample_name'], os.linesep,
				'>Consensus', os.linesep,
				dargs['iconsen'], '\t', dargs['sample_name'], os.linesep,
			])

with open('input_run.txt', 'a') as a:
	if dargs['igenes']:
		a.writelines(['>Genelist', os.linesep, dargs['igenes'], '\t', dargs['sample_name'], os.linesep])
	if dargs['iregion']:
		a.writelines(['>Region', os.linesep, dargs['iregion'], '\t', dargs['sample_name'], os.linesep])


# execute vprofiler Perl script
cmd = ' cat input_run.txt;   '

## required input
cmd += ' /usr/bin/perl %s/bin/vprofiler.pl -i input_run.txt -o Vpro ' % (cfg.customToolsDir, )

## optional input
if dargs['nt']:
	cmd += ' -nt '

if dargs['codon']:
	cmd += ' -codon -allheat '

if dargs['haplo']:
	cmd += ' -haplo '
	#cmd += ' -haplosep %s ' % (dargs['haplosep'])
	cmd += ' -haploseq %s ' % (dargs['haploseq'])
	if dargs['cq']:
		cmd += ' -consseq %s ' %  (dargs['consseq'])
		cmd += ' -conshaplo %s ' % (dargs['conshaplo'])

cmd += ' -noendvariant %s ' % (dargs['noendvariant'])

## output to dev/null since the output makes Galaxy break due to unicode string handling
cmd += ' &> /dev/null; '

#cmd = cmd + '/usr/bin/perl /home/gbotha/galaxy-dist/tools/custom_tools/bin/vprofiler.pl -i input_run.txt -o Vpro -nt -noendvariant=10 '
#cmd = cmd + '/usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/vprofiler.pl -i input_run.txt -o Vpro -nt -noendvariant=10 -codon -haplo -haploseq &> log;   '

## Zip the two output folders
cmd = cmd + ' zip -r Vpro_Heatmap_All.zip Vpro_%s_Heatmap_All &> /dev/null;  ' % (dargs['sample_name'])
cmd = cmd + ' zip -r Vpro_haplotypes.zip Vpro_haplotypes  &> /dev/null;  '

#'vprofiler.pl -i input.txt -o Vpro -noendvariant=10 -nt -codon -haplo -haploseq '
print cmd
subprocess.check_call(cmd, shell=True)

if os.path.exists('Vpro_' + dargs['sample_name'] + '_ntfreq.txt'):
	shutil.copyfile('Vpro_' + dargs['sample_name'] + '_ntfreq.txt', 'Vpro_ntfreq.txt')
if os.path.exists('Vpro_' + dargs['sample_name'] + '_codonfreq.xls'):
	shutil.copyfile('Vpro_' + dargs['sample_name'] + '_codonfreq.xls', 'Vpro_codonfreq.xls')
if os.path.exists('Vpro_' + dargs['sample_name'] + '_codon_details.txt'):
	shutil.copyfile('Vpro_' + dargs['sample_name'] + '_codon_details.txt', 'Vpro_codon_details.txt')

# these arent necesary 
#shutil.copyfile(dargs['out_ntfreq'], 'XXX_T0_ntfreq.txt')
#shutil.copyfile(dargs['out_codonfeq'], 'XXX_YYY_codonfreq.xls')
#shutil.copyfile(dargs['out_codondetails'], 'XXX_YYY_codon_details.txt')
#shutil.copyfile(dargs[''], 'XXX_YYY_')
#shutil.copyfile(dargs[''], 'XXX_YYY_')




















