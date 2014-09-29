#!/usr/bin/python

'''

Python wrapper for rc454.pl file. The wrapper was required since the the Perl script only takes a "basename" output option and builds output files from there. I had difficulties getting Galaxy to locate the correct files

'''

import os
import subprocess
import argparse
import pipes

# Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for rc454.pl, ensure that reads have been mapped \
						using Mosaik and that a QLX file exists')

## required args
parser.add_argument('qlx',     		help = 'Reads in QLX format')
parser.add_argument('fa',      		help = 'Reads in FASTA format')
parser.add_argument('qual',    		help = 'Read quality file in QUAL format')
parser.add_argument('ref',     		help = 'Reference sequence in FASTA format')

## optional args
parser.add_argument('--fullSettings')
parser.add_argument('--minhomosize')
parser.add_argument('--nqvalue')
parser.add_argument('--nocafie')
#parser.add_argument('--nqsmainqual')
#parser.add_argument('--nqsareaqual')
parser.add_argument('--nqssize')
parser.add_argument('--details')
parser.add_argument('--slicesize')
parser.add_argument('--gap3window')
parser.add_argument('--askPrimers')
parser.add_argument('--primers')
parser.add_argument('--primerbuffer')
parser.add_argument('--minpctvar')
parser.add_argument('--minnbvar')
parser.add_argument('--askGenes')
parser.add_argument('--genelist')
parser.add_argument('--noorf')

args = parser.parse_args()
dargs =  vars(args)
print dargs


## Perl file execution

cmd = ' /usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/rc454/rc454.pl '

## required inputs
cmd += ' %s %s %s %s ' % (dargs['qlx'], dargs['fa'], dargs['qual'], dargs['ref'])
## output prefix
cmd += ' rc454_out '

## optional inputs
if dargs['fullSettings'] == '1':
	cmd += ' -minhomosize %s' % dargs['minhomosize']
	cmd += ' -nqvalue %s ' % dargs['nqvalue']
	if dargs['nocafie'] == '1':
		cmd += ' -nocafie '
	#cmd += ' -nqsmainqual %s ' % dargs['nqsmainqual']
	#cmd += ' -nqsareaqual %s ' % dargs['nqsareaqual']
	cmd += ' -nqssize %s' % dargs['nqssize']
	if dargs['details'] == '1':
		cmd += ' -details '
	cmd += ' -slicesize %s ' % dargs['slicesize']
	cmd += ' -gap3window %s ' % dargs['gap3window']
	if dargs['askPrimers'] == '1':
		cmd += ' -primers %s ' % dargs['primers']
		cmd += ' -primerbuffer %s ' % dargs['primerbuffer']
	cmd += '-minpctvar %s ' % dargs['minpctvar']
	cmd += ' -minnbvar %s ' % dargs['minnbvar']
	if dargs['askGenes'] == '1':
		 cmd += ' -genelist %s ' % dargs['genelist']
	if dargs['noorf'] == '1':
		cmd += ' -noorf '  

## settings required for BROAD pipeline

cmd += '-bam'  

## zip if details are requested

if dargs['details'] == '1':
        cmd = cmd + '; zip -r rc454_out_details.zip rc454_out*  &> /dev/null;  '

#cmd = 'perl ~/galaxy-dist/tools/custom_tools/bin/rc454.pl %s %s %s %s  -bam ' % (pipes.quote(dargs['qlx']), pipes.quote(dargs['fa']), pipes.quote(dargs['qual']), pipes.quote(dargs['ref']))

print cmd

subprocess.check_call(cmd, shell=True)

#os.rename('XXX_cleaned.fasta', dargs['out_fasta'])
#os.rename('XXX_cleaned.qual', dargs['out_qual'])
#os.rename('XXX_final.sam', dargs['out_sam'])
#os.rename('XXX_final.bam', dargs['out_bam'])
#if os.path.exists(os.path.join(os.getcwd(), 'XXX_final_unmappedIDs.txt')):
#        os.rename('XXX_final_unmappedIDs.txt', dargs['out_txt'])
#if os.path.exists(os.path.join(os.getcwd(), 'XXX_final.qlx')):
#        os.rename('XXX_final.qlx', dargs['out_qlx'])
