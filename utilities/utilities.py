import os
import random
from itertools import islice

def ensure_dir_exists(x):
	if not os.path.exists(x):
		os.makedirs(x)

def write_cleaned_to_file(wpath, order, seqs):
    with open(wpath, 'w') as w:
        for seq in order:
            w.writelines(['>' + seq, os.linesep, seqs[seq], os.linesep])

def clean_seqs_fixed_lines(ifile, fastq=False):
	'''
		I'm trying a new approach here, where data is streamed line by line to be processed, instead of the entire file first being "slurped" up into memory.
	'''
	d		= {}
	order   = []
	if fastq:
		quals       = {}
		n=4
	else:
		n = 2
	with open(ifile, 'r') as f:
		last_position = f.tell()
		lines_gen = islice(f, n)
		lines = [line.strip() for line in lines_gen]
		while lines != []:
			if lines[0] == '':
				f.seek(last_position)
				f.readline()
			else:
				sname = lines[0][1:].strip().strip(os.linesep).strip('\r').strip()
				if sname in order:
					while sname in order:
						sname = sname + '_added-random_' + str(random.randint(0,10000))
				seq = lines[1].strip().strip(os.linesep).strip('\r').strip().upper()
				order.append(sname)
				d[sname] = seq
				if fastq:
					qual = lines[3].strip().strip(os.linesep).strip('\r').strip().upper()
					quals[sname] = qual
			last_position = f.tell()
			lines_gen = islice(f, n)
			lines = [line.strip() for line in lines_gen]
	
	if fastq:
		return order, d, quals
	else:
		return order, d


def clean_seqs(ifile, fastq=False):
	d		= {}
	q		= {}
	order	= []
	global sname
	if fastq:
		seq_mark = '@'
	else:
		seq_mark = '>'
	with open(ifile, 'r') as f:
		lines = f.readlines()
	#print lines
		for l in range(len(lines)):
			if lines[l].startswith(seq_mark):
				sname = lines[l][1:].strip().strip(os.linesep).strip('\r').strip()
				order.append(sname)
				tmp_lines = lines[l+1:]
				c = 0
				seq = ''
				try:
					while (not tmp_lines[c].startswith(seq_mark)) and (not tmp_lines[c] == os.linesep) and (not tmp_lines[c].startswith('+')):
						seq = seq + tmp_lines[c].strip().strip(os.linesep).strip('\r').strip().upper()
						c +=1
					#logger.debug('line section: %i \ntmp_lines: %s \ntmp_lines section: %i\n seq: %s ' % (l, tmp_lines, c, seq))
						if c == len(tmp_lines):
							break
				except IndexError:
					print tmp_lines, c, lines[l+1:l+5], sname

				if not sname in d.keys():
					d[sname] = seq.strip().strip(os.linesep).strip('\r').strip()
				else:
					x = order.index(sname)
					while sname in d.keys():
						sname = sname + '_added-random_' + str(random.randint(0,10000))
					order.remove(order[x])
					order.insert(x, sname)
					d[sname] = seq.strip().strip(os.linesep).strip('\r').strip()
			elif lines[l].startswith('+'):
				tmp_lines = lines[l+1:]
				c = 0
				seq = ''
				while (not tmp_lines[c].startswith(seq_mark)) and (not tmp_lines[c] == os.linesep) and (not tmp_lines[c].startswith('+')):
					seq = seq + tmp_lines[c].strip().strip(os.linesep).strip('\r').strip().upper()
					c +=1
					#logger.debug('line section: %i \ntmp_lines: %s \ntmp_lines section: %i\n seq: %s ' % (l, tmp_lines, c, seq))
					if c == len(tmp_lines):
						break
				sname = order[-1]
				if not sname in q.keys():
					q[sname] = seq.strip().strip(os.linesep).strip('\r').strip()
	if fastq:
		return order, d, q
	else:
		return order, d

def create_fasta(nams, seqs, wpath, s=None, e=None):
	with open(wpath, 'w') as w:
		for nam in nams:
			if s:
				w.writelines(['>', nam, '\n', seqs[nam][s:e], '\n'])
			else:
				w.writelines(['>', nam, '\n', seqs[nam], '\n'])

def create_fastq(nams, seqs, quals, wpath, s=None, e=None):
	with open(wpath, 'w') as w:
		for nam in nams:
			if s:
				w.writelines(['@', nam, '\n', seqs[nam][s:e], '\n', '+\n', quals[nam][s:e], '\n'])
			else:
				w.writelines(['@', nam, '\n', seqs[nam], '\n', '+\n', quals[nam], '\n'])
