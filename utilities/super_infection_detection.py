#!/usr/bin/python

'''

Author: Gordon

Date: 15 Oct 2013

The purpose of this script is to aid in detection super-infection detection. This is accomplished as follows:

1. Profile Align reads from T1 to panel and T2 to panel
2. Compare similarity between the (aligned) T2 reads and each of the panel sequences. Similarity is calculated in specific window sizes, which can vary.
2. (a) if any of the T2 sequences are sufficiently similar to the others, skip the process for those to speed up process. This "similarity threshold" must be settable by user
3. If similarity between T1 and T2 is LESS THAN the similarity between T2 and any of the panel sequences, we have a potential super-I eventi
3. (a) This is the start of the breakpoint. Move [jump_size] bases and repeat steps 1-3 until T1 seq no evidence of superI presents itself, which is the end of that breakpoint
4. Plot for each breakpoint window that was discovered. If T1 and T2 sequence is separated by more than one panel sequence (i.e. epidemiologically seperated), a super-I might be present.

Inputs:

1. Reference Panel Alignment (new or from history)
2. Reads from 2 time points, for any number of patients
3. Window Size (setting)
4. Similarity Threshold (setting)

'''
## TODO Find a way to throw out sequences ater trees have been created:
	## sort aoi's by size and iterate through them
	## use dendography to calc distances, compare distances to pick up broken trees
	## if no broken tree, try next aoi
	
## TODO OUTPUT
	## folder for each pair
		## tree file
		## output log highlighting findings
		## graph of similarities
		## fasta of alignment 

import os
import sys
import subprocess
import shutil
import numpy
import datetime
import cPickle as pickle
import utilities as util
import logging

import dendropy
from dendropy import treecalc

#x = os.path.join(dargs['od'], '002_logfile')
#logging.basicConfig(format= '%(asctime)-15s >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO, filename=x )

import maps
sim = maps.symmetrimalize(maps.RIP_lanl)

import argparse
parser = argparse.ArgumentParser(description = 'Detect super-Infectio events in reads from 2 different time points')

## Required input
parser.add_argument('--t1',  help='aligned reads from T1')
parser.add_argument('--t2',  help='aligned reads from T2')
parser.add_argument('--p',  help='panel reads to profile align to')
parser.add_argument('--s', help='similarity threshold to skip calcualtions from T2 that are too alike')
parser.add_argument('--w', help='window size of comparisons')
#parser.add_argument('--w', nargs='+', type=int, help='window size of comparisons')
#parser.add_argument('--am', default='clust', \
#		help='alignment method to use: \n clust -> clustalw, \n mos -> mosaik, \n ')
parser.add_argument('--j', help='size of window shift/jump for each sonsecutive similarity check')
parser.add_argument('--results', help='input results from file instead of calculating')
parser.add_argument('--stats', help='input stats from file instead of computing')
parser.add_argument('--od', help='name of the output directory where files will be written')
parser.add_argument('--draw', action="store_true", help='set whether trees must be draw')

args = parser.parse_args()
dargs =  vars(args)
print dargs

util.ensure_dir_exists(dargs['od'])
x = os.path.join(dargs['od'], '002_logfile')
logging.basicConfig(format= '%(asctime)-15s >> %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO,  filename=x )


##########################################################################################################

def similarity_check(s1, s2):
	''' 
Takes two sequences (or sequence sections) and computes a similarity score for them, based on the 'sim' matrix as imported \
from maps.py
Inputs:
sequences to be compared [STRING]
Returns:
similarity score [FLOAT]
	'''
        sim_score = 0
	logging.debug('			len s1: %i, len s2: %i' % (len(s1), len(s2)))
        if len(s1) != len(s2):
                logging.error("Length of sequences are not equal, make sure you've aligned your sequences. I have truncated to shortest seq")
		min_len = min(len(s1), len(s2))
        else:
		min_len = len(s1) 
        for i in range(min_len):
                try:
			#print float(sim[s1[i]][s2[i]])
                        sim_score += float(sim[s1[i]][s2[i]]) 
                except KeyError:
                        logging.debug('Could not calculate score on base "%s" in S1 and "%s" in S2 i \
                                - skipping position and moving forward' % (s1[i], s2[i]))
                        continue
	logging.debug("			sim_score : %f, len(s1) : %i, 100*(sim_score/min_len) : %f, dargs['s'] %f " % \
			(sim_score , len(s1), 100*(sim_score/min_len), float(dargs['s'])))
        
	return float("{0:.2f}".format(100*(sim_score/min_len)))


def profile_align(p, s):
	'''
Does a profile alignment for two sets of sequence. Each set should be aligned already. 
Inputs:
profile 1 an profile 2 [FILENAMES]
Returns:
Profiled sequences [DICT]
	'''
	logging.info('Performing Profile Alignment...')
        #cmd = ' /opt/exp_soft/clustalw-2.1/clustalw2 -profile1=%s -profile2=%s -output=fasta -outfile=profile.fasta -quiet' % (p, s) 
	cmd = ' /opt/exp_soft/muscle3.8/muscle3.8.31_i86linux64 -profile -in1 %s -in2 %s -out profile.fasta ' % (p, s)
	subprocess.check_call(cmd, shell=True)
        dummy, seqs_profiled = util.clean_seqs('profile.fasta')
        
        return seqs_profiled


def filter_similar(seqnams, seqs, seq_ref):
	'''
Compare all other sequences in a set to seq_ref. Sequences whose similarity is above the given threshold are filtered out of the \
seq_names list, which is used whenever looping through lists 
Inputs: 
sequence names in order [LIST]
sequences, stored under seq name [DICT]
name of sequence being compared to all others [STRING]
Returns:
Filtered version of sequence names [LIST]
	''' 
	logging.info('Filtering out sequences that are above theshold similarity (%s%%) for %s' % (dargs['s'], seq_ref))
	ss = list(seqnams)
	for i in range(len(seqnams)):
		s1 = seqs[seqnams[i]]
		logging.debug(' Comparing %s ...' % seqnams[i])
		j = i+1
		while j < len(seqnams):
			s2 = seqs[seqnams[j]]
			logging.debug('  with %s ... ' % seqnams[j])
			min_len = min(len(s1), len(s2))
			if similarity_check(s1, s2) > float(dargs['s']):
				if seqnams[j] in ss: 
					logging.info('   Sequences %s and %s are above similarity threshold, removing %s ' % \
					(seqnams[i], seqnams[j], seqnams[j]))
					ss.remove(seqnams[j])
			j +=1

	logging.info('Done...\n\n\n')	
	
	return ss


def compare_similarities(ws, results = {}):
	'''
For each sequence pair (one from time point 1 and one from time point 2), it compares the similarity between the sequences. \
If this similarity is lower than any of the similarity scores between the T2 sequence and the panel seqs, flag the area as \
an area of interest - which could potentially be a superinfection event.
There are a few intricacies in terms of extending an existing AOI and AOI's near the end of a sequence, but these steps are \
explained in the code.
Inputs:
None, all variables from global. Loops through all T1, T2 and panel sequences, for a given window size and jump size
Returns:
Results for each sequence pair where at least one AOI was found. Sim Scores for each T2 vs panel seq and T2 vs T1 are stored, \
as well as the AOI windows [DICT]
	'''
	logging.info('Comparing similarity scores for sequences to create areas of interest\n\n ')
	#results 		= {}
	ws			= int(ws)
	win_size		= ws
	js			= int(dargs['j'])	 
	for j, s2 in enumerate(t2_order_filtered):
		for i, s1 in enumerate(t1_order_filtered):
			number_steps 		= numpy.ceil( (len(p_t1_t2_seqs_profiled[s2])-ws) / float(js) ) +1 
			results[(s2, s1)] 	= {"scores" : {}, "aoi" : []}
			in_aoi 			= False
			end 			= False
			x 			= -js
			t2_vs_panel 		= numpy.zeros((number_steps, len(p_order)))
			t2_vs_t1    		= numpy.zeros(number_steps)
			while not end:
				x += js
				c = x/js
				if x >= len(p_t1_t2_seqs_profiled[s2])-ws:
					end = True
					ws  = len(p_t1_t2_seqs_profiled[s2])-x 
				mini = 1000000000
				for k, p in enumerate(p_order):
					s = t2_vs_panel[c, k] = similarity_check(p_t1_t2_seqs_profiled[s2][x:x+ws], \
						 p_t1_t2_seqs_profiled[p][x:x+ws])
					mini = min(s, mini)

				logging.debug('Pair: (%s, %s), Mini_panel : %f, x: %i, js: %i ' % (s2, s1, mini, x, js))
				t2_vs_t1[c] = t2t1 =  similarity_check(p_t1_t2_seqs_profiled[s2][x:x+ws], \
							p_t1_t2_seqs_profiled[s1][x:x+ws])
				logging.debug('Score (%s, %s) : %f ' % (s2, s1, t2t1))

				## last window at end of sequence
				if end:
					logging.debug(' Reached end of sequence.... x : %i, len(s1) : %i, ws : %i ' % (x, len(s1), ws))
					## if in AOI - end the AOI at the end of the sequence
					if in_aoi:
						logging.debug('		Ending current AOI...')
						finish = x + ws
						results[(s2, s1)]["aoi"].append((start, finish))
						#print results[(s2, s1)]["aoi"]
					## if not in AOI but a new one found, last AOI is set to window size
					elif t2t1 < mini:
						logging.debug('		Found new AOI in final window of sequence...')	
						start = x
						finish = x + ws
						results[(s2, s1)]["aoi"].append((start, finish))
						#print results[(s2, s1)]["aoi"]
					else:
						logging.debug('		Not currently in AOI and no new AOI found...')
				## intermediate window, minimum found
				elif t2t1 < mini:
					## start new AOI
					logging.debug(' SIMILARITY: a1 vs a2 is less similar than aribitrary panel similarity...')
					if not in_aoi:
						logging.debug('		Starting NEW AOI... ')
						in_aoi = True
						start = x
					## if already in AOI, continue with window
					else:
						logging.debug('		Already in AOI, continuing window size')
				## intermediate window, no minimum found
				else:
					## nothing happens keep going
					logging.debug(' No similarity interests in this window...')
					if not in_aoi:
						logging.debug('		Also not currently in AOI so moving on...')
					## end the current AOI and reset AOI marker
					else:
						logging.debug('		Ending current AOI')
						in_aoi = False
						finish = x + ws
						results[(s2, s1)]["aoi"].append((start, finish))
						#print results[(s2, s1)]["aoi"]
			t2_vs_panel_pickle_name = os.path.join(dargs['od'], '%s_vs_panel-jsize-%s_simthresh-%s_winsize-%s '.strip() % \
               				(s2, dargs['j'], dargs['s'], win_size))
			logging.debug('Pickling t2_vs_panel to file: %s' % (t2_vs_panel_pickle_name))
			dumpFile(t2_vs_panel, t2_vs_panel_pickle_name)
			results[(s2, s1)]["scores"][(s2, 'panel')] = t2_vs_panel_pickle_name	
			t2_vs_t1_pickle_name =  os.path.join(dargs['od'], '%s_vs_%s-jsize-%s_simthresh-%s_winsize-%s '.strip() % \
                                        (s2, s1, dargs['j'], dargs['s'], win_size))	
			logging.debug('Pickling t2_vs_t1 to file: %s' % (t2_vs_t1_pickle_name))
			dumpFile(t2_vs_t1, t2_vs_t1_pickle_name)
			results[(s2, s1)]["scores"][(s2, s1)] = t2_vs_t1_pickle_name
		if results[(s2, s1)]["aoi"] == []:
			logging.info('No AOI found for sequence pair: (%s, %s) - trashing...' % (s2, s1))
			results.pop((s2, s1))
		else:
			logging.info('%i AOI found for sequence pair: (%s, %s) ...' % ( len(results[(s2, s1)]["aoi"]), s2, s1))

	logging.info('Done...\n\n\n')
	
	return results

	
def compute_stats(res, ws):
	'''
Computes the following statistics for a given set:
1. Average similarity between sequences in a set for a given window size
2. Number / % of sequences with at least one AOI
3. Average number of AOI per sequence
4. Average AOI length 
	'''
	logging.info('Computing stats for data sets')
	stats		= {}
        ws      	= int(ws)
	win_size	= ws
        js      	= int(dargs['j'])
	## T1 average similarities
	logging.info(' T1 Average Similarity... (based on first 30 unfiltered sequences): ')
	try:
		tmp 		= min(len(t1_order_filtered), 30)
		t1_vs_t1_list 	= [] 
                for i in range(tmp-1):
                        for j in range(i+1, tmp):
				end = False
				x = -js
				c = x/js
				while not end:
					x += js
					if x >= len(p_t1_t2_seqs_profiled[t1_order_filtered[i]])-ws:
						end = True
						ws  = len(p_t1_t2_seqs_profiled[t1_order_filtered[i]])-x
					z = similarity_check(p_t1_t2_seqs_profiled[t1_order_filtered[i]][x:x+ws], p_t1_t2_seqs_profiled[t1_order_filtered[j]][x:x+ws])
					t1_vs_t1_list.append(z)	
		t1_vs_t1 = numpy.array(t1_vs_t1_list)
		#print t1_vs_t1
		stats["t1_ave"] = numpy.mean(t1_vs_t1)
		logging.info('\t %.2f %% ' % stats["t1_ave"])
		#t1_t1_pickle_name = os.path.join(dargs['od'], '%s_jsize-%s_simthresh-%s_winsize-%s '.strip() % \
		#	('t1_vs_t1', dargs['j'], dargs['s'], win_size))
		#stats["t1_scores"] = t1_t1_pickle_name
		#dumpFile(t1_vs_t1, t1_t1_pickle_name)
		t1_vs_t1_list   = []
		t1_vs_t1	= []	
		
	except MemoryError:
		logging.info(' MemoryError ')

	## T2 average similarities
        logging.info(' T2 Average Similarity... (based on first 30 unfiltered sequences) ')
        try:
                tmp             = min(len(t2_order_filtered), 30)
                t2_vs_t2_list   = [] 
                for i in range(tmp-1):
                        for j in range(i+1, tmp):
                                end = False
                                x = -js
                                c = x/js
                                while not end:
                                        x += js
                                        if x >= len(p_t1_t2_seqs_profiled[t2_order_filtered[i]])-ws:
                                                end = True
                                                ws  = len(p_t1_t2_seqs_profiled[t2_order_filtered[i]])-x
                                        z = similarity_check(p_t1_t2_seqs_profiled[t2_order_filtered[i]][x:x+ws], p_t1_t2_seqs_profiled[t2_order_filtered[j]][x:x+ws])
                                        t2_vs_t2_list.append(z)  
                t2_vs_t2 = numpy.array(t2_vs_t2_list)
                stats["t2_ave"] = numpy.mean(t2_vs_t2)
		logging.info('\t %.2f %% ' % stats["t2_ave"])
                #t2_t2_pickle_name = os.path.join(dargs['od'], '%s_jsize-%s_simthresh-%s_winsize-%s '.strip() % \
                #        ('t2_vs_t2', dargs['j'], dargs['s'], win_size))
                #stats["t2_scores"] = t2_t2_pickle_name
                #dumpFile(t2_vs_t2, t2_t2_pickle_name) 
		t2_vs_t2_list   = []
		t2_vs_t2	= []

        except MemoryError:
                logging.info(' MemoryError ')

	## Panel Average Similarities
        logging.info(' P Average Similarity... (based on first 30 unfiltered sequences) ')
        try:
                tmp             = min(len(p_order), 30)
                p_vs_p_list   = [] 
                for i in range(tmp-1):
                        for j in range(i+1, tmp):
                                end = False
                                x = -js
                                c = x/js
                                while not end:
                                        x += js
                                        if x >= len(p_t1_t2_seqs_profiled[p_order[i]])-ws:
                                                end = True
                                                ws  = len(p_t1_t2_seqs_profiled[p_order[i]])-x
                                        z = similarity_check(p_t1_t2_seqs_profiled[p_order[i]][x:x+ws], p_t1_t2_seqs_profiled[p_order[j]][x:x+ws])
                                        p_vs_p_list.append(z)  
                p_vs_p = numpy.array(p_vs_p_list)
                stats["p_ave"] = numpy.mean(p_vs_p)
		logging.info('\t %.2f %% ' % stats["p_ave"])
                #p_p_pickle_name = os.path.join(dargs['od'], '%s_jsize-%s_simthresh-%s_winsize-%s '.strip() % \
                #        ('p_vs_p', dargs['j'], dargs['s'], win_size))
                #stats["p_scores"] = p_p_pickle_name
                #dumpFile(p_vs_p, p_p_pickle_name)
		p_vs_p_list 	= []
		p_vs_p		= []

        except MemoryError:
                logging.info(' MemoryError ')


	## Additional Stats
	logging.info(' Additional Stats... ')
	try:
		## Number and percentage of T2 sequences with AOI's
		pairs = [x for x in res.keys() if type(x) == tuple] 
		stats["seqs_with_aoi"] 		= list(set([ x[0] for x in pairs]))
		stats["num_seqs_with_aoi"]	= len(stats["seqs_with_aoi"])
		#print (len(stats["seqs_with_aoi"]),  len(t2_order))
		stats["perc_seqs_with_aoi"]	= 100* (len(stats["seqs_with_aoi"]) / float(len(t2_order_filtered)))
		logging.info(' %.2f %% of the T2 sequences had AOI' % stats["perc_seqs_with_aoi"])
		
		## Average number of AOI per sequence pair
		stats["total_num_AOI"]		= {}
		stats["total_num_AOI_list"]	= []
		for s1 in t1_order_filtered:
			for s2 in t2_order_filtered:
				if (s2, s1) in pairs:
					stats["total_num_AOI"][(s2, s1)] = len(res[(s2, s1)]["aoi"])
					stats["total_num_AOI_list"].append(len(res[(s2, s1)]["aoi"]))
		stats["ave_aoi_per_seq_pair"]	= numpy.average(stats["total_num_AOI_list"])
		logging.info(' Average number of AOI per sequence pair (of those who had at least one AOI): %f ', \
			  stats["ave_aoi_per_seq_pair"])

		## Average AOI length
		stats["aoi_length_per_seq_pair_list"] = []
		for keypair in pairs:
			for aoi in res[keypair]["aoi"]:
				#print aoi
				#print keypair
				#print res[keypair]["aoi"]
				len_aoi = int(aoi[1]) - int(aoi[0])
				stats["aoi_length_per_seq_pair_list"].append(len_aoi)
		stats["ave_aoi_length"] = numpy.average(stats["aoi_length_per_seq_pair_list"])
		logging.info(' Average length of AOI: %.2f ', stats["ave_aoi_length"])
	except KeyError:
		logging.info(' Error...')

	return stats

def compare_trees(results):
	'''
	loops through the AOI's for each seq pair and computes tress for the ares. If s2 and s1 are epidemiologically seperated in the tree, the AOI is kept in results, else it is popped
	'''
	logging.info(' Combing through trees to find epidemiological splits:...')
	#key = results.keys()[0]
	#print key
	#stuffs = [key]
	#print 'fff', results.keys()
	z = [d for d in results.keys() if type(d)==tuple]
	#print 'zzz', z
	for x  in z:
		logging.info('\t checking pair (%s, %s)' % (x[0], x[1]))
		odir = os.path.join(dargs['od'], '005_trees_of_AOI_before_evaluation_for_epi_splits', x[0], x[1])
		util.ensure_dir_exists(odir)
		windows = sorted(results[x]['aoi'], key= lambda item : (item[1]-item[0]), reverse=True)
		epi_split_found = False
		counter = -1
		while not epi_split_found:
			counter +=1
			window = windows[counter]
			logging.info('\t\t window... (%s, %s)' % (window[0], window[1]))
			## indent the block below back to this level if you want all windows and remove above loop
			#print window 		
			s = window[0]
			e = window[1]
			#print p_t1_t2_seqs_profiled[x[1]][window[0]: window[1]]
			#print p_t1_t2_seqs_profiled[x[0]][window[0]: window[1]]
			wname = '_'.join([x[0], 'vs', x[1], str(window[0]) + '-' + str(window[1]), '.fasta'])
			wpath = os.path.join(odir, wname)
			nams = [x[0]] + [x[1]] + p_order
			#print nams
			util.create_fasta( nams, p_t1_t2_seqs_profiled, wpath, s, e )
			fpath = wpath
			wname = '_'.join([x[0], 'vs', x[1], str(window[0]) + '-' + str(window[1]), '.tree'])
			wpath = os.path.join(odir, wname)
			cmd = 'FastTreeMP -nt < %s > %s -quiet' % (fpath, wpath)
			subprocess.check_call(cmd,shell=True)	
			
			# cycle through trees, calculate distances for pair seqs and all seqs paired with them.
			fpath = wpath
			tree = dendropy.Tree.get_from_path(fpath, 'newick')
			pdm = treecalc.PatristicDistanceMatrix(tree)
			for i, node in enumerate(tree.taxon_set):
				if node.label.replace(' ', '_') == x[0]:
					for node2 in tree.taxon_set[i+1:]:
						mini = 1
						pair_dist = 10000
						if node2.label.replace(' ', '_') == x[1]:
							pair_dist = pdm(node, node2)
						else:
							mini = min(pdm(node, node2), mini)
							mini_seq = node2.label.replace(' ', '_')
					if mini == pair_dist:
						if window == windows[-1]:
							logging.info(' Pair: ( %s, %s) is not EPI split in final window: (%i, %i) - removing from results' % \
								(x[0], x[1], window[0], window[1])) 
							results.pop(x)
						else:
							logging.info(' Pair: ( %s, %s) is not EPI split in window: (%i, %i) - moving to next window' % \
								(x[0], x[1], window[0], window[1]))

					else:
						logging.info(' An EPI split has been detected for pair: ( %s, %s) in window: (%i, %i) - plot will be drawn. ' % \
							(x[0], x[1], window[0], window[1]))
						logging.info('\t The pair was split by sequence: %s ' % (mini_seq)) 
						epi_split_found = True
						copy_dir = os.path.join(dargs['od'], 'final_results', x[0], x[1])
						util.ensure_dir_exists(copy_dir)
						shutil.copy(os.path.join(odir, '_'.join([x[0], 'vs', x[1], \
							str(window[0]) + '-' + str(window[1]), '.fasta'])), \
							copy_dir)
						shutil.copy(os.path.join(odir, '_'.join([x[0], 'vs', x[1], \
                                                        str(window[0]) + '-' + str(window[1]), '.tree'])), \
                                                        copy_dir)

	return results

def dumpFile(x, wpath):
	'''
Dumps 'x' into a file in wpath by python pickling it. This is to enable future runs of the same data without \
redoing the alignments and similarity checks
	'''
	#dump_file_name =  os.path.join(dargs['od'], '%s_%s_jsize-%s_simthresh-%s_winsize-%s '.strip() % \
	#	(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S_'), prefix, dargs['j'], dargs['s'], ws)) 
	with open(wpath, 'wb') as w:
		pickle.dump(x, w)	


###################################################################################################################3
## MAIN ##

## create output dir if it doesn't exist
if not dargs['od']:
	sys.exit('You need to specify the --od option for the output directory')
else:
	out_align = os.path.join(dargs['od'], '003_cleaned_alignments')
	util.ensure_dir_exists(out_align)
	out_pickles = os.path.join(dargs['od'], '004_pickles')
	util.ensure_dir_exists(out_pickles)

if not dargs['results']:

	results = {}

	## read seqs
	logging.info('Reading Sequences...')
	t1_order, t1_seqs       = util.clean_seqs(dargs['t1'])
	t1_cleaned              = os.path.join(out_align, \
					dargs['t1'].split('.')[0] + '_CLEANED.' + dargs['t1'].split('.')[1])      
	util.write_cleaned_to_file(t1_cleaned, t1_order, t1_seqs)
	
	t2_order, t2_seqs       = util.clean_seqs(dargs['t2'])
	t2_cleaned              = os.path.join(out_align, \
					dargs['t2'].split('.')[0] + '_CLEANED.' + dargs['t2'].split('.')[1])
	util.write_cleaned_to_file(t2_cleaned, t2_order, t2_seqs)
	
	p_order, p_seqs 	= util.clean_seqs(dargs['p'])
	p_cleaned               = os.path.join(out_align, \
					dargs['p'].split('.')[0] + '_CLEANED.' + dargs['p'].split('.')[1])
	util.write_cleaned_to_file(p_cleaned, p_order, p_seqs)
	logging.info('Done...\n\n\n')
	#print p_order[0:2], p_seqs[p_seqs.keys()[0]]


	## profile align T1 and T2 reads to Panel
	p_and_t1_seqs_profiled 	= profile_align(p_cleaned, t1_cleaned)
	p_and_t1_cleaned	= os.path.join(out_align, 'panel_t1_profile.fasta')
	util.write_cleaned_to_file(p_and_t1_cleaned, t1_order + p_order, p_and_t1_seqs_profiled) 
	logging.info('Created profile alignment for T1')
	
	p_t1_t2_seqs_profiled = profile_align(p_and_t1_cleaned, t2_cleaned)
	p_t1_t2_cleaned        = os.path.join(out_align, 'panel_t1_t2_profile.fasta')
	util.write_cleaned_to_file(p_t1_t2_cleaned, t2_order + t1_order + p_order, p_t1_t2_seqs_profiled)
	logging.info('Created profile alignment for T2') #: %s, Alignment: \n %s ' % ( t2_order[0], t2_seqs_profiled[t2_order[0]]))


	## check for similarity between input sequences to remove sequences that are too similar to the rest.
	#print sim
	t1_order_filtered = filter_similar(t1_order, p_t1_t2_seqs_profiled, 'T1')
	t2_order_filtered = filter_similar(t2_order, p_t1_t2_seqs_profiled, 'T2')

	results['t1_order_filtered'] = t1_order_filtered

	results['t2_order_filtered'] = t2_order_filtered
	
	results['p_order'] = p_order


ws = dargs['w']
## Compare similarity scores for window of size dargs['w'] for T2 vs T1 seqs, moving window by jump size (dargs['j']) 
## Store each score in the results dictionary under the key (t2,t1) and the sub-key (t2,px)
## Store windows where similarities differ under the key (t2,t1) and the sub-key "windows", which is a list of (start, end) tuples
## If no windows exist for a pair, delete the keys from the results dictionary (use pop) 
if dargs['results']:
	try:
		logging.info('Importing results from file %s ' % dargs['results'])
		with open(dargs['results'], 'rb') as f:
			results = pickle.load(f)
		t1_order_filtered 	= results['t1_order_filtered']
		t2_order_filtered 	= results['t2_order_filtered']
		p_order			= results['p_order']
		dummy, p_t1_t2_seqs_profiled = util.clean_seqs(os.path.join(out_align, 'panel_t1_t2_profile.fasta') )
	except:
		logging.error('Could not retrieve results from the given file. Typo? Is the file in this working DIR?')
else:
	results = compare_similarities(ws, results)
	## dump to file
	dump_file_name =  os.path.join(out_pickles, 'results_winsize-%s '.strip() %  ws)
	dumpFile(results, dump_file_name)
#print results.keys()
#print results[('CAP281_313_044_F9', 'CAP281-1')]["aoi"]

## Compute statistics & dump to file
if dargs['stats']:
	try:
		logging.info('Importing stats from file %s ' % dargs['stats'])
		with open(dargs['stats'], 'rb') as f:
			stats = pickle.load(f)
	except:
		logging.error('Could not retrieve stats from the given file. Typo? Is the file in this working DIR?')
else:
	stats = compute_stats(results, ws)
	dump_file_name = os.path.join(out_pickles, 'stats_winsize-%s '.strip() %  ws)
	dumpFile(stats, dump_file_name)

with open(os.path.join(dargs['od'], '001_summary'), 'w') as w:
	w.writelines(["\n\nAverage similarities for T1, T2 and Panel:\n\n", str(stats["t1_ave"]), \
		'\t', str(stats["t2_ave"]), '\t', str(stats["p_ave"])])
	w.writelines(["\n\n %.2f %% of the T2 sequences had AOI\n\n" % stats["perc_seqs_with_aoi"]])
	w.writelines(['\n\nAverage number of AOI per sequence pair (of those who had at least one AOI): %f\n\n ' \
		% stats["ave_aoi_per_seq_pair"]])
	w.writelines(['\n\n Average length of AOI: %.2f \n\n'% stats["ave_aoi_length"]])
	w.writelines(["\n\n T2 Sequences with epidemiologically split AOI's:\n\n"])
	w.writelines([str(x) + '\n' for x in stats["seqs_with_aoi"]])
	
#print stats.keys()
print stats["t1_ave"], stats["t2_ave"], stats["p_ave"]
print stats["seqs_with_aoi"], stats["num_seqs_with_aoi"], stats["perc_seqs_with_aoi"]
#print stats["total_num_AOI_list"]
print stats["ave_aoi_per_seq_pair"]
#print stats["ave_aoi_length"]

## Build trees for significant results
## create input file
#print results.keys()				
#print results[results.keys()[0]]


#with open(os.path.join(dargs

if dargs['draw']:
	results = compare_trees(results)			

allowed_keys 		= [d for d in results.keys() if type(d)==tuple]
sorted_allowed_keys 	= sorted(allowed_keys, key = lambda item : (item[0])) 
with open(os.path.join(dargs['od'], '001_summary'), 'a') as w:
	w.writelines(["\n\nThe following (T2, T1) pairs had epidemiologically split AOI's in them (final_results folder): \n\n"])
	w.writelines(["(" + str(x[0]) + ", " + str(x[1]) + ')\n' for x in sorted_allowed_keys])
#z = sorted(z,  key= lambda item : (item[0]))
for i, s2 in enumerate(t2_order_filtered):
	for j, s1 in enumerate(t1_order_filtered):
		if (s2, s1) in allowed_keys: 	
			# save data as csv
			t2_vs_t1 = pickle.load(open(results[(s2, s1)]["scores"][(s2, s1)], 'rb'))
			t2_vs_panel = pickle.load(open(results[(s2, s1)]["scores"][(s2, 'panel')], 'rb'))
			#print dir(t2_vs_t1)
			#print t2_vs_t1.size
			#print t2_vs_panel.size
			odir = os.path.join(dargs['od'], 'final_results', s2, s1)
			util.ensure_dir_exists(odir)
			wname = "%s_vs_%s.csv" % (s2, s1)
			wpath = os.path.join(os.path.join(odir, wname))
			#t2_vs_t1.tofile(open(wpath, 'w'), sep=",", format="%s")
			with open(wpath, 'w') as w:
				w.write('T2_seq, Window, %s [T1],  '% s1)
				for seq in p_order[:-1]:
					w.write('%s [P], ' % seq)
				w.write('%s [P]\n ' % p_order[-1])
				for c in range(t2_vs_panel.shape[0]):
					w.write('%s, ' % s2)
					w.write('%s, ' % str(c+1))
					w.write('%s, ' % t2_vs_t1[c])
					for v in range(t2_vs_panel.shape[1]-1):
						w.write('%s, ' % t2_vs_panel[c,v])
					w.write('%s\n' % t2_vs_panel[c,-1])
				
			
			#wname = "%s_vs_%s.csv" % (s2, 'panel')
			#wpath = os.path.join(os.path.join(odir, wname))
			#t2_vs_panel.tofile(open(wpath, 'w'), sep=",", format="%s")
		
		
			# if any other pair is smaller than seq pair, flag the result as meaningful
				# save data as csv
				# create R plot script and save
				# plot the thing and output it

			# if not, remove from results and log it
		
			

## TODO use the UPGMA algo for tree-building (faster)
## TODO the PIM option outputs the percentage identity matrix for faster comparisons??



















