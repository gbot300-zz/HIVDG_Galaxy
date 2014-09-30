#!/usr/bin/python

'''

Author: Gordon

Date 13 Dec 2013

This tool sorts sequence reads by a specific sequence polymer in the sequence (usually denoted as the PID). It can be used to sort multiplexed reads by sample / patient / experiment etc.

Essentially, the tool takes as input a pre-PID sequence (can be degenerate) and a post-PID sequence (can also be degenerate) and searches for these polymers in each read. The PID is then identified (if it is present) and the read is grouped with the other reads from that PID. 

A consensus sequences is then generated from the reads in each PID bin.

One is able to specify whether a sample of 1000 bins and their reads needs to be output to a zipped folder for QC purposes
'''

## TODO

## 1. Presently, consensus sequences can be aligned as an option. This option should be removed to ensure more modularity.
## 2. The script runs very slow due to high memory reqs - can be made much more memory efficient -> http://effbot.org/zone/wide-finder.htm

## Below is an example of what the PID and PRE and POST strings will look like in a read.

#V3 reverse primer sequence (in rev comp ie: as it should appear in the final merged file)

#CAGGAGGGGAYCTAGAArTTACAACnnnnnnnnnCTGAGCGTGTGGCAAGGC

#gene specific primer		PID		   universal primer site

#V1 reverse primer sequence (in rev comp ie: as it should appear in the final merged file)

#CACYTKAGAATCGCATAACCAGCTGGnnnnnnnnnCTGAGCGTGTGGCAAGGC

#gene specific primer		PID		   universal primer site

import os
import subprocess
import sys

import utilities as util

import maps

import random

#from fuzzywuzzy import fuzz, process
import re

import logging
#logging = logging.getLogger('Logger')
#logging.setLevel('INFO')

import argparse
parser = argparse.ArgumentParser(description = 'Cycles through a set of merged Illumina reads and identifies Primer IDs. New PIDs are compared to an existing set and if mismatched by only one base, the sequence is added to the existing PID')

## Required input
parser.add_argument('--fa', help='Merged reads in fasta format')
parser.add_argument('--fq', help='Merged reads in fastq format')
parser.add_argument('--rfa', help='Reference sequence in fasta format, if left blank, first sequence in reads file will be used as reference')
parser.add_argument('--rfq', help='Reference sequence in fastq format, if left blank, first sequence in reads file will be used as reference')
parser.add_argument('--pre', help='Pre PID sequence - degeneracy allowed')
parser.add_argument('--PID', help='PID sequence - most of the time should be "N"s only')
parser.add_argument('--post', help='Post Primer sequence - degeneracy not allowed')
parser.add_argument('--sample_name', help='Name of sample to be prefixed to consensus alignments')
parser.add_argument('--consensus', action="store_true", help='Do you want to find a consensus sequence for the sequences in each bin?')
parser.add_argument('--cut', action="store_true", help='Do you want the PID and trailing bases to be cut? - PID will be moved to sequence name regardless')
parser.add_argument('--lines', action="store_true", help='Is the fasta/fastq file in standard format, each sequence on a line (as opposed to sequence split over multiple lines')
parser.add_argument('--dump', action="store_true", help='For intermediate steps, dump sequences fopr each bin into a file, which can be rerun through this tool')
parser.add_argument('--sample_alignments', action="store_true", help="Store the first 100 alignments for the bins for QC purposes")
parser.add_argument('--min_upids', help="minimum number of unidentical representative sequences required to create a consensus sequence for a PID")
parser.add_argument('--min_ipids', help="minimum number of identical representative sequences required to create a consensus sequence for a PID")
parser.add_argument('--wtd_nopids', help="What to do with sequences in which no PIDS were detected") 
parser.add_argument('--wtd_opids', help="What to do with sequences that have less than the required number of PIDs for consensus")
parser.add_argument('--pid_summary', action="store_true",  help="Do you want to output a CSV format file with the number of sequences to each PID?")

args = parser.parse_args()
dargs =  vars(args)
#print dargs

logFilename = "log-" + str(random.randint(0,10000))
logging.basicConfig(format= '%(asctime)-15s || %(levelname)s  >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', \
filename=logFilename, filemode='w', level=logging.INFO)

base_replacer = {
		'A' : '[A|M|R|W|D|H|V|N]',
		'C' : '[C|M|S|Y|B|H|V|N]',
		'G' : '[G|R|S|K|B|D|V|N]',
		'T' : '[T|W|Y|K|B|D|H|N]'
		}


## FUNCTION DEFINITIONS

def create_consensus(aligned_seqs, max_len):
''' 
This is an adaptation of a function written by Ram from SANBI. See "Make_Consensus.py" in the Seq2Res code
RETURNS:
	Consensus sequence for a set of aligned reads
'''
		#max_len=len(consensus_list[0][1])
		position_collection_box={}
		consensus_sequence=""
		for position in range(max_len):
			## first build the frequencies for each base at each position
			for seq in aligned_seqs.keys():
				readsequence = aligned_seqs[seq]
				if len(readsequence)>position:
					base=readsequence[position]
					if base in position_collection_box.keys():
						position_collection_box[base]+=1
					else:
						position_collection_box[base] =1
				else:
					break
			## secondly, determine the maximum frequency base. if there are ties create lists
			max_freq = -1
			for key in	position_collection_box.keys():
				if position_collection_box[key] > max_freq:
					base_to_add = key
					max_freq = position_collection_box[key]
				elif position_collection_box[key] == max_freq:
					if isinstance(base_to_add, list):
						base_to_add.append(key)
					else:
						base_to_add = [base_to_add]
						base_to_add.append(key)
			## if a list of candidates exists, look up the ambiguity character
			if isinstance(base_to_add, list):
				## hack here until we figure out what to do with gaps. 
				base_to_add = sorted([x.upper() for x in base_to_add])
				if (len(base_to_add) == 2) and ('-' in base_to_add):
					base_to_add = '-'
				else:
					bta_list = base_to_add[:]
					for x in bta_list:
						if x not in ['A', 'C', 'G', 'T']:
							base_to_add.remove(x)
					for key in maps.translator.keys():
						base_to_add.sort()
						if maps.translator[key] == base_to_add:
							base_to_add = key

			#print base_to_add
			try:
				consensus_sequence+=base_to_add
			except:
				base_to_add = 'N'
			position_collection_box.clear()
		return consensus_sequence

def build(fixed, vary, build_len):
'''
This function builds all possible patterns of length [build_len], from a sequence [vary], which may contain degenerate bases - hence the call to maps.translator. This creates all possible variants of a given PRE / POST input string w.r.t. degenaracy
RETURNS:
	all possible variants of a given length for a degenerate sequence string
'''
	for i in range(len(vary)):
		fx = fixed[:]
		for fix in fx:
			for j in maps.translator[vary[i]]:
				fixed.append(fix+j)
		fixed = [x for x in fixed if len(x) == i+1]
	builds = [x for x in fixed if len(x) == build_len]
				
	return builds

def build_patterns(pres, pid, posts):
'''
Given a set of PRE patterns and POST patterns, this function tries to find all combinations of PRE and POST in a sequence in order to identify the PID. The base_replacer functionality enables one to match ACGT in the PRE / POST string to degenerate bases in the read.
RETURNS:
	all possible PRE-PID-POST pattern combinations for the given PRE, PID and POST patterns
'''
	patterns = []
	for post in posts:
		for pre in pres:
			pattern = '%s(%s)%s' % (pre, '.'*len(PID), post)
			logging.debug('pattern added: %s' % (pattern))
			for base in ['A', 'C', 'G', 'T']:
				pattern = pattern.replace(base, base_replacer[base])
			patterns.append(pattern)
	return patterns

def id_pids(seqs):
'''
Iterates through the reads and tries to identify PIDs. If a new one is found, the PID is added as a key to the unique_pids dictionary object and the read is stored under the PID. If the PID already exists, the read is added to the current reads under that PID. If no PID is identified the read is recorded in the no_pids list.
RETURNS: 
unique_pids dict object with upids as keys and all read names of reads in that upid bin.
no_pids list of read names which contain no PID
'''
		for sname in seqs.keys():
				#print sname
				## first step is to identify the primer seq
				matched = False
				c = 0
				logging.debug('\nsname %s\nseq %s\npre %s\npost %s' % (sname, seqs[sname], pre, post))
				for pattern in patterns:
						c+=1
						pid_seq = re.compile(pattern)
						if pid_seq.search(seqs[sname]):
								pid = pid_seq.search(seqs[sname]).group(1)
								logging.debug('MATCHED:   ::: %s' % (pid))
								matched = True
								if not pid in unique_pids.keys():
										unique_pids[pid] = [sname]
								else:
										unique_pids[pid].append(sname)
								break
				if not matched:
						logging.debug('NO MATCH')
						no_pids.append(sname)
		return unique_pids, no_pids

def get_max_len(order, aligned_seqs):
'''
obtain the max length of all the reads for a given PID. This value is used during consensus genration
RETURNS: max_len of reads integer value
'''
	max_len = 0
	for seq in order:
		#print len(aligned_seqs[seq]), max_len
		if len(aligned_seqs[seq]) > max_len:
			max_len = len(aligned_seqs[seq])
	return max_len


def align_and_consensus(con_keys, upids):
'''
For each unique PID, [upid], a file is created containing all the reads for that upid and the reads are aligned using mafft.
A consensus sequence for each PID bin is then generated from the aligned reads in each bin.
If the sample alignments option has been enabled, the first 1000 bins' reads alignments are kept and output as a zipped folder for QC purposes
RETURNS: a disctionary object with consensus names as keys and consensus sequences as values 
'''
		con_seqs = {}
		randy = str(random.randint(0,1000000))
		if dargs['sample_alignments']:
			util.ensure_dir_exists('sample_alignments')
			samples = 0
		else:
			samples = 1000
		for upid in con_keys:
			tmpnam		= 'tmp' + randy 
			util.create_fasta(upids[upid], seqs, tmpnam) #(read_names, reads, filename) essentially. 
			mafftnam	= 'mafft_out_msa' + randy ## TODO remember when you change this ramics to change the dump aligns for nopids and opids too
			subprocess.check_call("mafft --auto --quiet %s	> %s;" % (tmpnam, mafftnam), shell=True)
			order, aligned_seqs = util.clean_seqs(mafftnam)
			max_len = get_max_len(order, aligned_seqs)
			con_seq		= create_consensus(aligned_seqs, max_len)
			conseq_nam	= dargs['sample_name'] + '_' + upid + '_' + str(len(upids[upid]))
			con_seqs[conseq_nam] = con_seq
			if samples < 1000:
				subprocess.check_call('cp %s sample_alignments/%s.fasta;' % (mafftnam, conseq_nam), shell=True)
				samples +=1
		#print con_seqs.keys()[-1:], con_seqs[con_seqs.keys()[1]], con_seqs[con_seqs.keys()[2]]
		return con_seqs

def create_fasta_custom(nams, seqs, wpath, s=None, e=None):
'''
Writes out sequences to a a file in fasta format. One is able to write only parts of the sequences using the 's' and 'e' options
'''
		with open(wpath, 'w') as w:
				for nam in nams:
						if s:
							w.writelines(['>', nam, '\n', seqs[nam][s:e], '\n'])
						else:
							w.writelines(['>', nam, '\n', seqs[nam], '\n'])


def dump_bins_to_file(upids, seqs, quals=[]):
'''
outputs the reads for a number of unique PIDS into a file. This is used for QC purposes with the [sample_alignments] argument
'''
	# create_fastq(nams, seqs, quals, wpath, s=None, e=None):
	bindir = "binned" + str(random.randint(0,100000))
	os.mkdir(bindir)
	for upid in upids:
		if dargs['fq']:
			fname	=	upid + ".fastq"
			fpath	=	os.path.join(bindir, fname)
			util.create_fastq(upids[upid], seqs, quals, fpath)
		else:
			fname	=	upid + ".fasta"
			fpath	=	os.path.join(bindir, fname)
			util.create_fasta(upids[upid], seqs, fpath)
	subprocess.check_call(" zip -r binned.zip %s &> /dev/null; " % (bindir), shell=True)
	

#######################################################################################3
if __name__ == "__main__":
	logging.info("Reading in data...")

	## used to specify whether the input file is in the old Sanger sequence format that splits reads over multiple new lines
	if dargs['lines']:
		read_seqs = util.clean_seqs_fixed_lines
	else:
		read_seqs = util.clean_seqs

	## read in sequences, from fasta or fastq
	if dargs['fa']:
		logging.info("Fasta file detected")
		snames, seqs		= read_seqs(os.path.join(os.getcwd(), dargs['fa']), fastq=False) 
	else:
		logging.info("FastQ file detected")
		snames, seqs, quals		= read_seqs(os.path.join(os.getcwd(), dargs['fq']), fastq=True)

	#print snames
	#print seqs
	#print quals

	## read in reference sequence, from fasta or fastq
	if dargs['rfa'] :
		ref, rseq,		= util.clean_seqs(os.path.join(os.getcwd(), dargs['rfa']), fastq=False)
	elif dargs['rfq'] :
		ref, rseq, rqual	= util.clean_seqs(os.path.join(os.getcwd(), dargs['rfq']), fastq=True)
	else:
		sys.exit('Must enter a reference sequence')

	logging.info('Done reading data')
	## build the possible combinations of the pre-primer sequence
	pre		= dargs['pre'].upper()
	logging.info("Building PRE patterns...")
	pre_builds	= build([''], pre, len(pre))
	logging.info("done")

	## build the possible combos of the post-primer sequence. this isnt stricly necessary as these aren't degenrate but might as well
	post		= dargs['post'].upper()
	logging.info("Building POST patterns...")
	post_builds		= build([''], post, len(post))
	logging.info("done\n\n")

	PID		= dargs['PID'].upper()

	## build RE patterns to match for in each seq
	patterns = build_patterns(pre_builds, PID, post_builds)

	## identify unique PIDS, identify seqs without PIDs
	unique_pids				= {}
	no_pids					= []
	unique_pids, no_pids	= id_pids(seqs)

	## minimum number of unidentical and identical sequences required to create consensus
	if dargs['consensus']:
		min_upids				= int(dargs['min_upids'])
		min_ipids				= int(dargs['min_ipids'])
	else:
		min_upids = min_ipids = 1

	## log some info
	num_seqs				= len(seqs.keys())
	num_nopids				= len(no_pids)
	num_upids				= max(len(unique_pids.keys()),1)

	logging.info('Total number of unique PID detected: %i' % (num_upids))

	#ones					= [x for x in unique_pids.keys() if len(unique_pids[x]) == 1]
	#num_ones				= len(ones)
	#logging.info('Number of PIDs represented by only 1 sequence: %i (%.2f%%)' % (num_ones, 100* float(num_ones)/num_upids))

	# seqs[unique_pids[x][0]] != seqs[unique_pids[x][1]]
	#len(set([seqs[unique_pids[x][i] for i in len(unique_pids[x]]]))
	if dargs['consensus']:
		unidentical_pids = []
		for x in unique_pids.keys():
			if len(unique_pids[x]) >= min_upids:
				tmp = [seqs[unique_pids[x][i]] for i in range(len(unique_pids[x]))]
				if len(set(tmp)) > 1:
					unidentical_pids.append(x)
		num_unidentical_pids	= len(unidentical_pids) 
		logging.info('Number of PIDs represented by %s or more (unidentical) sequences: %i (%.2f%%)' % (str(min_upids), num_unidentical_pids, 100* float(num_unidentical_pids)/num_upids))

		identical_pids = []
		for x in unique_pids.keys():
			if len(unique_pids[x]) >= min_ipids:
				tmp = [seqs[unique_pids[x][i]] for i in range(len(unique_pids[x]))]
				if len(set(tmp)) == 1:
					identical_pids.append(x)
		num_identical_pids		= len(identical_pids)
		logging.info('Number of PIDs represented by %s or more identical sequences: %i (%.2f%%)' % (str(min_ipids), num_identical_pids, 100* float(num_identical_pids)/num_upids))

		z						= identical_pids + unidentical_pids + no_pids 
		other_pids				= [x for x in unique_pids.keys() if not x in z]
		other_pids_seqs			= []
		for x in other_pids:
			other_pids_seqs = other_pids_seqs + unique_pids[x]

	logging.info('Total number of sequences: %i ' % (num_seqs))

	perc					= 100* float(num_nopids) / num_seqs
	logging.info('Number of sequences in which no PID were detected: %i (%.2f%%) ' % (num_nopids, perc))

	if dargs['consensus']:
		numseqs_otherpids		= len(other_pids_seqs) 
		perc2					= 100* float(numseqs_otherpids) / num_seqs
		logging.info('Number of sequences in which a PID was detected but representation criteria was not met: %i (%.2f%%) ' % (numseqs_otherpids, perc2)) 

		pids					= unidentical_pids + identical_pids 
		numseqs_pids			= sum([len(unique_pids[x]) for x in pids])
		perc3					= 100* float(numseqs_pids) / num_seqs
		logging.info('Number of sequences in which a PID was detected and representation criteria was met: %i (%.2f%%)\n\n ' % (numseqs_pids, perc3))
		## what to do with  noPID sequences and sequences below thresholds above
		#if dargs['wtd_nopids'] == "


	## if this is an intermediate step, store all the sequences for a UID in a file and zip all these files together
	## otherwise align seqs for each primer & create consensus from that alignment (if specified)

	if dargs['dump']:
		logging.info('Creating intermediate bin directory')
		if dargs['fq']:
			dump_bins_to_file(unique_pids, seqs, quals)
		else:
			dump_bins_to_file(unique_pids, seqs)
		logging.info('done...')
	elif dargs['consensus']:
		logging.info('Creating consensus files...')
		consensus_seqs	= {}
		consensus_keys	= unidentical_pids + identical_pids
		num_conkeys		= len(consensus_keys)
		consensus_seqs	= align_and_consensus(consensus_keys, unique_pids)
		
		logging.info('Consensus file being written in FASTA format')
		create_fasta_custom(sorted(consensus_seqs.keys()), consensus_seqs, 'consensus_sequences.fasta')
		logging.info('done...')
		if dargs['sample_alignments']:
			logging.info('zipping sample alignments')
			subprocess.check_call( 'zip -r sample_alignments.zip sample_alignments &> /dev/null;', shell=True)
			logging.info('done...')

	logging.info('Handling reads with no PID')
	if dargs['wtd_nopids'] == "dump":
		if dargs['fa']:
			util.create_fasta(no_pids, seqs, "nopids.fasta")
		elif dargs['fq']:
			util.create_fastq(no_pids, seqs, quals, "nopids.fastq")
	elif dargs['wtd_nopids'] == "align":
		if dargs['fa']:
			util.create_fasta(no_pids, seqs, "nopids.fasta")
			subprocess.check_call("mafft --auto --quiet nopids.fasta  > nopids_aligned.fasta", shell=True)
		elif dargs['fq']:
			util.create_fastq(no_pids, seqs, quals, "nopids.fastq")
			subprocess.check_call("mafft --auto --quiet nopids.fastq  > nopids_aligned.fastq", shell=True)
	logging.info('done...')

	logging.info('Handling reads with PID"s present, but that did not meet representative requirements')
	if dargs['wtd_opids'] == "dump":
		if dargs['fa']:
			util.create_fasta(other_pids_seqs, seqs, "other_pids.fasta")
		elif dargs['fq']:
			util.create_fastq(other_pids_seqs, seqs, quals, "other_pids.fastq")
	elif dargs['wtd_opids'] == "align":
		if dargs['fa']:
			util.create_fasta(other_pids_seqs, seqs, "other_pids.fasta")
			subprocess.check_call("mafft --auto --quiet other_pids.fasta > other_pids_aligned.fasta", shell=True)
		elif dargs['fq']:
			util.create_fastq(other_pids_seqs, seqs, quals, "other_pids.fastq")
			subprocess.check_call("mafft --auto --quiet other_pids.fastq > other_pids_aligned.fastq", shell=True)
	logging.info('done...')

	if dargs['pid_summary']:
		logging.info('Creating PID Summary file')
		sx = sorted(unique_pids.items(), key=lambda z : len(z[1]), reverse=True)
		with open('pid_summary.csv', 'w') as w:
			w.writelines(['PID, number of sequences, sequence names\n'])
			for item in sx:
				pid			= item[0]
				sequences	= item[1]
				pid_seq_nams = ",".join([nam for nam in sequences])
				w.writelines([pid, ',', str(len(sequences)), ',', pid_seq_nams , '\n'])
		logging.info('done')


	logging.info('moving LOG file')
	subprocess.check_call('mv %s log.txt;' % (logFilename), shell=True)
	logging.info('done...')








