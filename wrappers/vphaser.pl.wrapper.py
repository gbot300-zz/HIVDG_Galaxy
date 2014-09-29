#!/usr/bin/python

'''

Python wrapper for vphaser.pl file. The wrapper was required since the the Perl script only takes a "basename" output option and builds output files from there. I had difficulties getting Galaxy to locate the correct files

'''

import shutil
import subprocess
import argparse

# Define usage options using argparse library
parser = argparse.ArgumentParser(description = 'Wrapper for vphaser.pl, in preparation for nucleotide, \
						codon and haplotype identification with V-Profiler')
## required arguments
parser.add_argument('--i',      		help = 'QLX format input file from RC454 aligner')

## optional arguments
parser.add_argument('--fullSettings')
parser.add_argument('--model')
parser.add_argument('--qqplot')
parser.add_argument('--offset')
parser.add_argument('--l')

args = parser.parse_args()
dargs =  vars(args)
print dargs

# Perl file execution
# I had to add this as Galaxy was keep the file as a '.dat' file, which made the vphaser script break
cmd = ' /usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/vphaser.pl '

## required args
## this hack is required since the vphaser.pl tool wants a file named "something.qlx"
shutil.copyfile(dargs['i'], 'zzz.qlx') 
cmd += ' -i zzz.qlx '
cmd += ' -o vphaser_out '

## optional args
if dargs['fullSettings'] == '1':
	cmd += ' -model %s ' % (dargs['model'])
	cmd += ' -offset %s ' % (dargs['offset'])
	cmd += ' -l %s ' %(dargs['l'])
	if dargs['qqplot'] == '1':
		cmd += ' -q '

print cmd

cmd += ' &> /dev/null '
 
subprocess.check_call(cmd, shell=True)

#shutil.copyfile(dargs['i'], 'zzz.qlx')
#cmd = '/usr/bin/perl /opt/exp_soft/galaxy-tools/custom_tools/bin/vphaser.pl -i zzz.qlx -o xyz'




