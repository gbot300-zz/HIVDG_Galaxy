#!/usr/bin/python
#get and compare AA frequencies bwtn 2 alignments using fisher's exact test
#test for residue association where pvalue from above <0.05
#to run:  /opt/rocks/bin/python aafreqs.py aafreqsSensitive.txt aafreqsResistant.txt aafreqPVALUES aafreqCOUNTS

import sys
from rpy import *
import os
infile1,infile2, outfile, logfile = sys.argv[-4], sys.argv[-3],sys.argv[-2], sys.argv[-1]			#tab_del aa alignment files

def aafrequency(alignment1, alignment2):			#for each alignment
	f1 = open(alignment1,"r")
	lines1 = f1.readlines()
	sub1 = []				#list of file1 seqs
	allresults = {}
	
	for line1 in lines1:
		lin = line1.strip().split("\t")
		ID = lin[0].strip()
		seq = lin[-1].strip()
		sub1.append(seq)
	
	f2 = open(alignment2,"r")
        lines2 = f2.readlines()
        sub2 = []                               #list of file2 seqs
	
        for line2 in lines2:
                lin2 = line2.strip().split("\t")
                ID2 = lin2[0].strip()
                seq2 = lin2[-1].strip()
                sub2.append(seq2)

	for i in range(len(sub1[0])):			#for each pstn in alignment
		s = str(i+1)			#actual site
		allaa = []			#list of all aa at this pstn from both C and B
		for sq in sub1:
			if sq[i] not in allaa:	#add all res to allaa list
				allaa.append(sq[i])
		for sq in sub2:
			if sq[i] not in allaa:
				allaa.append(sq[i])
	#NOW GET THE aa  in each alignment at this site
		afreqs = {}	#aa:[count1,count2]	
		aa1 = []	#list of actual freqs in 1 for stats test
		aa2 = []	#list of actual freqs in 2 for stats test
		for a in allaa:	#get freq of each res 
			afreqs[a] = [0,0]
			for ssq in sub1:
				if ssq[i] == a:
					afreqs[a][0] = afreqs[a][0]+1
			for ssq in sub2:
				if ssq[i] == a:
					afreqs[a][-1] = afreqs[a][-1]+1
	#	print i+1, afreqs		
		
		for r in afreqs:
			aa1.append(afreqs[r][0])
			aa2.append(afreqs[r][-1])
	#	print aab, aac
		allresults[str(i+1)] = [afreqs,aa1,aa2]
	f1.close()
	f2.close()
	allresults = sorted(allresults.items())
	return allresults

data = aafrequency(infile1, infile2)
outf = open(outfile,"w")
outf.write("Site\tPvalue\tPvalue(BonferroniCorrected)\n")
outf2 = open(logfile,"w")
outf2.write("Site\tResidue\tCountINfile1\tCountINfile2\n")

numtests = 0			#get num of tests to be done
for nt in data:
	if len(nt[-1][1])>1:	#sites with one and identical residue btwn alingments not tested
		numtests = numtests+1

for sit in range(len(data)):
	site = data[sit][0]
#	print data[sit]
	sited = data[sit][-1]
#	print sited
	for aa in sited[0]:		#get res count per site into outf2
#		print aa, sited[0][aa][0]
		outf2.write(site+"\t"+aa+"\t"+str(sited[0][aa][0])+"\t"+str(sited[0][aa][-1])+"\n")

	if len(sited[1])>1:		#do fisher's exact test on aa counts
		r_fisher_file = os.environ.get("R_FISHER_FILE")
		r.source(r_fisher_file)
		#print sited[1],sited[-1]
		fisher = r.my_fisher(sited[1],sited[-1])
#		if fisher['p.value'] < 0.045:
#			print site, data[site][1], data[site][-1],"fisher pvalue = : ", fisher['p.value'], fisher['p.value']*len(data)
		outf.write(site+"\t"+str(fisher['p.value'])+"\t"+str(fisher['p.value']*numtests)+"\n")
		if fisher['p.value'] < 0.045:
			outf.write("AA\tpvalue\tBonferroniCorrected\tOddsRatio\n")
			for sa in sited[0]:
				sa1 = [sited[0][sa][0],int(sum(sited[1])-sited[0][sa][0])]	#count aa in file1, count rest aa
				sa2 = [sited[0][sa][-1],int(sum(sited[-1])-sited[0][sa][-1])]
#				print sa1, sa2
				fisherassoc = r.my_fisher(sa1,sa2)
				outf.write(sa+"\t"+str(fisherassoc['p.value'])+"\t"+str(fisherassoc['p.value']*len(sited[0]))+"\t"+str(fisherassoc['estimate']['odds ratio'])+"\n")
	else:
		outf.write(site+"\t"+"similar_single_residue"+"\n")
outf.close()
outf2.close()
