#!/usr/bin/python

'''

Python wrapper for vphaser2 file. 

'''

import shutil
import subprocess
import argparse
import glob
import cfg as cfg

# Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for vphaser2, in preparation for nucleotide, \
		                        codon and haplotype identification with V-Profiler')
## required arguments
parser.add_argument('--i',              help = 'BAM format')

## optional arguments
parser.add_argument('--fullSettings', help = 'Full settings switch')
parser.add_argument('--e',	help = 'pile-up and phasing switch')
parser.add_argument('--w', help = 'alignment_window size')
parser.add_argument('--ig', help='number of bases to ignore on each end')
parser.add_argument('--delta', help='constrain PE distance by delta x fragsize_variation (auto measured by program)')
parser.add_argument('--ps', help='default 30; percentage of reads to sample to get stats.')
parser.add_argument('--dt', help='default 1; 1: dinucleotide for err prob measure; 0: not')
parser.add_argument('--cy', help='default 1; 1: read cycle for err calibr; 0: not')
parser.add_argument('--mp', help=' default 1; 1: mate-pair for err calibr; 0: not')
parser.add_argument('--qual', help='default 20; quantile of qual for err calibr')
parser.add_argument('--a', help='default 0.05; significance value for stat test')
parser.add_argument('--ref', help='The reference file required for the vprofiler format creator')

args = parser.parse_args()
dargs =  vars(args)
print dargs

cmd = ' variant_caller '

cmd += ' -i %s -o . ' % ( dargs['i'])

if dargs['fullSettings']:
	cmd += ' --e %s ' % dargs['e']
	cmd += ' --w %s ' % dargs['w']
	cmd += ' --ig %s  ' % dargs['ig']
	cmd += ' --delta %s  ' % dargs['delta']
	cmd += ' --ps %s  ' % dargs['ps']
	cmd += ' --dt %s ' % dargs['dt']
	cmd += ' --cy %s ' % dargs['cy']
	cmd += ' --mp %s ' % dargs['mp']
	cmd += ' --qual %s ' % dargs['qual']
	cmd += ' --a %s ' % dargs['a']

print cmd

cmd += ' &> /dev/null '

subprocess.check_call(cmd, shell=True)

fdrvar	= [x for x in glob.glob('*.fdr.var.txt')][0]
covplot = [x for x in glob.glob('*.covplot.R')][0]

shutil.copy(fdrvar, 'vphaserout.txt')
shutil.copy(covplot, 'covplot.R')

subprocess.check_call( 'R CMD BATCH covplot.R', shell=True)

cmd = '/usr/bin/perl %s/bin/vph2vprf_format.pl -vph2 vphaserout.txt -ref %s -o vpro.txt' % (cfg.customToolsDir, dargs['ref'])
subprocess.check_call( cmd , shell=True)

