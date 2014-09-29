#!/usr/bin/python

'''

Python wrapper for runMosaik2.pl file. The wrapper was required since the the Perl script only takes a "basename" output option and builds output files from there. I had difficulties getting Galaxy to locate the correct files

'''

import os
import subprocess
import argparse

# Define usage options using argparse library
# compulsory args
parser = argparse.ArgumentParser(description = 'Wrapper for runMosaik2.pl, which is a precursor for using the RC454 tood')
parser.add_argument('--st',	help = 'Sequencing tech used')
parser.add_argument('--fa', 	help = 'Reads in FASTA format')
parser.add_argument('--qual', 	help = 'Read quality file in qual format')
parser.add_argument('--fa2',      help = 'Reads in FASTA format (second pairs if paired reads)')
parser.add_argument('--qual2',    help = 'Read quality file in qual format (second pairs if paired reads')
parser.add_argument('--fq',      help = 'Reads in FASTQ format')
parser.add_argument('--fq2',      help = 'Reads in FASTQ format (second pairs if paired reads')
parser.add_argument('--ref', 	help = 'Reference sequence in FASTA format')
# default args (only add these if default settings are selected) 
#parser.add_argument('paramillu',	help = 'Recommended parameter settings for Illumina data')
#parser.add_argument('param454',	help = 'Recommended parameter settings for 454 data')
#parser.add_argument('qlx', 		help = 'Required for running RC454')
# advanced args
parser.add_argument('--hs', help ='')
parser.add_argument('--act', help ='')
parser.add_argument('--mmp', help ='')
parser.add_argument('--minp', help ='')
parser.add_argument('--ms', help ='')
parser.add_argument('--mms', help ='')
parser.add_argument('--gop', help ='')
parser.add_argument('--hgop', help ='')
parser.add_argument('--gep', help ='')
parser.add_argument('--nqsmq', help ='')
parser.add_argument('--nqsaq', help ='')
parser.add_argument('--nqssize', help ='')
parser.add_argument('--bw', help ='')
parser.add_argument('--nqvalue', help ='')
parser.add_argument('--m', help ='')
parser.add_argument('--mfl', help ='')
#parser.add_argument('out_sam', help = 'Name of the QLX format output file')
#parser.add_argument('out_sam', help = 'Name of the SAM format output file')
#parser.add_argument('out_bam', help = 'Name of the BAM format output file')
#parser.add_argument('out_txt', help = 'Name of the TXT format output file')

args = parser.parse_args()
dargs =  vars(args)

#print dargs

# Perl file execution
#cmd = '/usr/bin/perl /home/gbotha/galaxy-dist/tools/custom_tools/bin/runMosaik2.pl -fa %s -qual %s -ref %s -o runMos_out -qlx -param454' % (dargs['fa'], dargs['qual'], dargs['ref'])
cmd = '/usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/runMosaik2.pl'

# seq tech
cmd += ' -st %s ' % dargs['st']

# input files
if dargs['fa']:
	if dargs['fa2']:
		cmd += ' -fa %s -fa2 %s -qual %s -qual2 %s ' % (dargs['fa'], dargs['fa2'], dargs['qual'], dargs['qual2'])
	else:
		cmd += ' -fa %s -qual %s ' % (dargs['fa'], dargs['qual'])
elif dargs['fq']:
	if dargs['fq2']:
		cmd += ' -fq %s -fq2 %s ' % (dargs['fq'], dargs['fq2'])
	else:
		cmd += ' -fq %s '  % (dargs['fq']) 

# ref seq
cmd += ' -ref %s ' % (dargs['ref'])

# output name
cmd += ' -o runMos_out -nqsmq 20 -nqsaq 15 ' # Hack here for the NQSMQ and NQSAQ options that need to be present when creating a QLX file

# default options
if not dargs['hs']:
	if dargs['st'] == '454':
		cmd += ' -param454 '
	elif dargs['st'] == 'illumina':
		cmd += ' -paramillu '

# advanced options
else:
	cmd += ' -hs %s' % dargs['hs']
	cmd += ' -act %s' % dargs['act']
	cmd += ' -mmp %s' % dargs['mmp']
	cmd += ' -minp %s' % dargs['minp']
	cmd += ' -ms %s' % dargs['ms']
	cmd += ' -mms %s' % dargs['mms']
	cmd += ' -gop %s' % dargs['gop']
	cmd += ' -hgop %s' % dargs['hgop']
	cmd += ' -gep %s' % dargs['gep']
	cmd += ' -nqsmq %s' % dargs['nqsmq']
	cmd += ' -nqsaq %s' % dargs['nqsaq']
	cmd += ' -nqssize %s' % dargs['nqssize']
	cmd += ' -bw %s' % dargs['bw']
	cmd += ' -nqvalue %s' % dargs['nqvalue']
	cmd += ' -m %s' % dargs['m']
	if dargs['mfl']:
		cmd += ' -mfl %s' % dargs['mfl']

# qlx required for RC454
cmd += ' -qlx '

#cmd = '/usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/runMosaik2.pl -fa %s -qual %s -ref %s -o runMos_out -qlx -param454' % (dargs['fa'], dargs['qual'], dargs['ref'])

print cmd
subprocess.check_call(cmd, shell=True)

#os.rename('.sam', dargs['out_sam'])
#os.rename('.bam', dargs['out_bam'])
#if os.path.exists(os.path.join(os.getcwd(), '.qlx')):
#	os.rename('.qlx', dargs['out_qlx'])
#if os.path.exists(os.path.join(os.getcwd(), 'XXX_unmappedIDs.txt')):
#	os.rename('_unmappedIDs.txt', dargs['out_txt'])




 
