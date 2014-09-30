#!/usr/bin/python

'''

Galaxy wrapper for the bin_by_seqstring.py tool. This tool search a input sequence file for a specific pattern in each sequence, in order to bin the sequences by Primer ID's, sample ID's etc.
The tool is set up in such a way that one can repetitively bin the sequences for consecutive runs.

'''

import os
import subprocess
import argparse
import cfg as cfg

## Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for bin_by_seqstring.py')

## required args
parser.add_argument('--fa', help='Merged reads in fasta format')
parser.add_argument('--fq', help='Merged reads in fastq format')
parser.add_argument('--rfa', help='Reference sequence in fasta format, if left blank, first sequence in reads file will be used as reference')
parser.add_argument('--rfq', help='Reference sequence in fastq format, if left blank, first sequence in reads file will be used as reference')
parser.add_argument('--pre', help='Pre PID sequence - degeneracy allowed')
parser.add_argument('--PID', help='PID sequence - most of the time should be "N"s only')
parser.add_argument('--post', help='Post Primer sequence - degeneracy not allowed')
parser.add_argument('--sample_name', help='sample name to be prefixed to consensus sequence names')

## optional args
parser.add_argument('--consensus', action="store_true", help='Do you want to find a consensus sequence for the sequences in each bin?')
parser.add_argument('--cut', action="store_true", help='Do you want the PID and trailing bases to be cut? - PID will be moved to sequence name regardless')
parser.add_argument('--lines', action="store_true", help='Is the fasta/fastq file in standard format, each sequence on a line (as opposed to sequence split over multiple lines')
parser.add_argument('--dump', action="store_true", help='Dump intermediate bins to a folder')
parser.add_argument('--sample_alignments', action="store_true", help="Store the first 100 alignments for the bins for QC purposes")
parser.add_argument('--min_upids', help="minimum number of unidentical representative sequences required to create a consensus sequence for a PID")
parser.add_argument('--min_ipids', help="minimum number of identical representative sequences required to create a consensus sequence for a PID")
parser.add_argument('--wtd_nopids', help="What to do with sequences in which no PIDS were detected")
parser.add_argument('--wtd_opids', help="What to do with sequences that have less than the required number of PIDs for consensus")
parser.add_argument('--pid_summary', action="store_true", help='Print out a CSV file summarising the number of sequences per PID and what they are')

args = parser.parse_args()
dargs =  vars(args)
print dargs

## tool execution
cmd = '/usr/bin/python %s/galaxy_confs/bin_by_seqstring.py' % (cfg.customToolsDir,) 

## required inputs
if dargs['fa']:
	cmd += ' --fa %s ' % (dargs['fa'])
elif dargs['fq']:
	cmd += ' --fq %s ' % (dargs['fq'])

if dargs['rfa']:
	cmd += ' --rfa %s ' % (dargs['rfa'])
elif dargs['rfq']:
	cmd += ' --rfq %s ' % (dargs['rfq'])

cmd += ' --pre %s --PID %s --post %s ' % (dargs['pre'], dargs['PID'], dargs['post'])

if dargs['dump']:
	cmd += ' --dump '
elif dargs['consensus']:
	cmd += ' --consensus '
	cmd += ' --sample_name %s' % (dargs['sample_name'])
	cmd += ' --min_upids %s ' % (dargs['min_upids'])
	cmd += ' --min_ipids %s ' % (dargs['min_ipids'])
	cmd += ' --wtd_opids %s ' % (dargs['wtd_opids'])

if dargs['cut']:
	cmd += ' --cut '

if dargs['lines']:
	cmd += ' --lines '

if dargs['sample_alignments']:
	cmd += ' --sample_alignments '

cmd += ' --wtd_nopids %s ' % (dargs['wtd_nopids'])

if dargs['pid_summary']:
	cmd += ' --pid_summary '

print cmd
subprocess.check_call(cmd, shell=True)


















