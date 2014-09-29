#!/usr/bin/python

'''

Python wrapper for vphaser.pl file. The wrapper was required since the the Perl script only takes a "basename" output option and builds output files from there. I had difficulties getting Galaxy to locate the correct files

'''

import shutil
import subprocess
import argparse
parser = argparse.ArgumentParser(description = 'Detect super-Infectio events in reads from 2 different time points')

parser.add_argument('--t1',  help='aligned reads from T1')
parser.add_argument('--t2',  help='aligned reads from T2')
parser.add_argument('--p',  help='panel reads to profile align to')
parser.add_argument('--s', help='similarity threshold to skip calcualtions from T2 that are too alike')
parser.add_argument('--w', help='window size of comparisons')
parser.add_argument('--j', help='size of window shift/jump for each sonsecutive similarity check')

args = parser.parse_args()
dargs =  vars(args)
print dargs



# Python file execution
#cmd = ' /usr/bin/python /home/gbotha/galaxy-dist/tools/custom_tools/bin/super_infection_detection.py '
cmd = ' /usr/bin/python /opt/exp_soft/galaxy-tools/custom_tools/bin/super_infection_detection.py '

## required args
## this hack is required since the vphaser.pl tool wants a file named "something.qlx"
cmd += ' --t1 %s ' % dargs['t1']
cmd += ' --t2 %s ' % dargs['t2']
cmd += ' --p %s ' % dargs['p']
cmd += ' --s %s ' % dargs['s']
cmd += ' --w %s ' % dargs['w']
cmd += ' --j %s ' % dargs['j']

cmd += ' --od . '
cmd += ' --draw '

print cmd

cmd += ' &> /dev/null ;'

## zip required files
cmd = cmd + ' zip -r 003_cleaned_alignments.zip 003_cleaned_alignments &> /dev/null; '
cmd = cmd + ' zip -r 005_trees_of_AOI_before_evaluation_for_epi_splits.zip 005_trees_of_AOI_before_evaluation_for_epi_splits &> /dev/null; '
cmd = cmd + ' zip -r final_results.zip final_results &> /dev/null; '
subprocess.check_call(cmd, shell=True)



