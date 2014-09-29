#!/usr/bin/python
'''

Author: Nobu

#search for presence of  LANL epitopes (>=70% identity, for specific HLAs) in first sequence of an Amino Acid tab_del file eg 1st sequence could be a vaccine_insert for example, then determines divergence (using HIV-Between distance matrix) of autologous sequences from the 1st sequence at the identified epitope sites for each HLA allele genotype in the infected host: Output is mean epitope distance for each patient: the results can be used to compare these distances between 2 sets of data eg between placebo and vaccinees, using a non-parametric equivalent of a t-test

#USER INPUT = "tab del amino acid file (1st sequence assumed to be the reference) & "tab del HLA data file, no title line" & "gene name, should be either Gag, Nef or Pol"

#e.g. TO RUN: python matchLANL.py lanl_seqtest.txt lanl_pHLA_all.txt Gag lanl_out1 lanl_out2 

#OUTPUT:- prints a list of PIDs, each PID with mean distance of allepitopes and another for optimal epitopes  
#OUTPUT; - file 2 has a list of HLAs each with a summary of predicted epitope start sites 

Edited by: Gordon Sep 2014
1. TODO Added Fasta input functionality
2. Wrote the transform_coords function to enable epitope searching in sequences that did not contain the entire HIV genome

'''

import sys, re, os

seqfile, hlafile, gene, outfile1, outfile2, = sys.argv[-5], sys.argv[-4], sys.argv[-3],sys.argv[-2], sys.argv[-1]
rel_start_pos = rel_end_pos = 0
#seqfile, hlafile, gene, outfile1, outfile2, rel_start_pos, rel_end_pos = sys.argv[-7], sys.argv[-6], sys.argv[-5],sys.argv[-4], sys.argv[-3], int(sys.argv[-2]), int(sys.argv[-1])		

mfile = open("/opt/exp_soft/galaxy-tools/custom_tools/meta_files/HIVbetween.txt","r")              #get matrix and AA codes
mlines = mfile.readlines()
aacodes = {}                            #residue, code A=0 etc
matrix = {}                             #aa:distscores for matrix row
for ml in range(len(mlines)):
        mline = mlines[ml].strip().split(",")
        aacodes[mline[0].strip()] = ml
        matrix[mline[0].strip()] = mline[1:]

def transform_coords(lines):
	'''
	Function added by Gordon to tranform the coords in the LANLallepitopes.txt file so that you can search for epitopes in your
	sequence, even if you only have a part of the HIV genome
	'''
	new_lines = []
	for line in lines:
		excl = False
		## extract the epitope positions from each line
		start, end = int(line.split('\t')[2]), int((line.split('\t')[3]))
		## tranform positions
		s = start - rel_start_pos
		e = end - rel_start_pos
		## if epitope is in HIV genome before where our sequence starts, exclude it
		if s < 0:
			excl = True
		## if epitope is in a region of the HIV genome after our sequences, exclude it
		if end > rel_end_pos:
			excl = True

		## replace line, if the epitope might be present
		if not excl:
			new_line = '\t'.join([line.split('\t')[0], line.split('\t')[1], str(s), str(e), line.split('\t')[4], line.split('\t')[5]])
			new_lines.append(new_line)

	#return new_lines
	return lines

def matchLANL(seqfile, hlafile, gene):

#get pid HLA data
	phla_list = []			#list all patient HLAs together
	allphla = {}                    #pid:[hla list]
	hlaf = open(hlafile,"r")	#get patient HLA data
	hlines = hlaf.readlines()
	for hline in hlines:
        	hlin = hline.strip().split("\t")
	        pd = hlin[0].strip()
        	phlas = hlin[1:]
		allphla[pd] = []
        	for phla in phlas:
	#		hhh = re.sub("\*", "", phla.strip())
			hhh = phla.strip()
	                if hhh not in allphla[pd]:
        	                allphla[pd].append(hhh)
                	else:
                        	pass
			if hhh not in phla_list:
				phla_list.append(hhh)
	hlaf.close()
#get seqdata, 1st sequence in alignment is reference/vaccine insert
	file = open(seqfile,"r")
	lines = file.readlines()
	vaccinen = lines[0].strip().split("\t")[0].strip()
	vaccines = re.sub("-","X",lines[0].strip().split("\t")[-1].strip())		#reference sequence
	data = {}		#seqdata
	for line in lines[1:]:
        	lin = line.strip().split("\t")
        	sname = lin[0].strip()
	        seq = re.sub("-","X",lin[-1].strip())
        	data[sname] = seq
	file.close()
#get LANL epitope list: HLA: [(startsite, epitopeseq),(),(),...,]
 
	lanfile = open("/opt/exp_soft/galaxy-tools/custom_tools/meta_files/LANLallepitopes.txt")			#LANLfile is either 'all' or 'alist' epitope file
	lanlines = lanfile.readlines()
	'''
    Gordon added the transform_coords function below to adjust the positional coordinates in the LANLallepitopes.txt file to enable 
	epitope searching if you dont have the entire HIV genome
	'''
	lanlines = transform_coords(lanlines)
	lanlhlas = {}				#HLA: [(startsite, epitopeseq),(,),(),....] for the named gene
	for lanline in lanlines:
		if gene in lanline:		#only lines for the named gene
			lanl = lanline.strip().split("\t")
			pep = lanl[0].strip()
			start = lanl[2].strip()
			hlass = lanl[6:]
			hlas = []
			for hla in hlass:
        			hl = hla.strip().split(",")
	                	for b in hl:
        	        		if b.strip() not in hlas:
                        #	print b.strip()
                	                	hlas.append(b.strip())
			for h in hlas:
				if h not in lanlhlas:
					lanlhlas[h] = [(start,pep)]
				else:
					if (start,pep) not in lanlhlas[h]:
						lanlhlas[h].append((start,pep))
					else:
						pass
 	lanfile.close()
#	print lanlhlas
#	print phla_list
	opt_lanlhlas = {}				#HLA: [(startsite, epitopeseq),(,),(),....] for the named gene
	optf = open("/opt/exp_soft/galaxy-tools/custom_tools/meta_files/LANLalistepitopes.txt")
	olines = optf.readlines()
	for oline in olines:
		if gene in oline:
			ohlas = []
			opep = oline.strip().split("\t")[0]
			ohlass = oline.strip().split("\t")[2:]
			for ohla in ohlass:
				ohl = ohla.strip().split(",")
				for o in ohl:
					if o.strip() not in ohlas:
						ohlas.append(o.strip())
			for oh in ohlas:
				if oh in lanlhlas:
					for pp in lanlhlas[oh]:
						if opep in pp:		#get opt peptide data
							if oh not in opt_lanlhlas:
								opt_lanlhlas[oh]=[pp]
							else:
								if pp not in opt_lanlhlas[oh]:
									opt_lanlhlas[oh].append(pp)
								else:
									pass
	optf.close()
#	print opt_lanlhlas 
				
#search for epitopes maching >= 70% on the ref seq for all PID hla alleles
	all_matches = {}			#hla:[[PID list][peptide_startsite list]]
	opt_matches = {}			#for optimal HLAs alone
	Pall_matches = {}			#PID:[list opt predicted epitope sites (unique)]
	Pall_matchesAB = {}			#FOR HLA-A and B alleles alone
	Popt_matches = {}
	Popt_matchesAB = {}
	for ph in allphla:			#for each individual
		for hh in allphla[ph]:
			def match70(lanllist, matches, Pmatches, AB):	#lanlhlas, all_matches Pallmatches, Pall_matchesAB OR opt_lanlhlas, opt_matches, Popt_matches, Popt_matchesAB
				if hh in lanllist:
		#			print hh
					for phd in lanllist[hh]:
						st = phd[0]		#start_site
						lp = phd[-1]		#LANL epitope sequence
						m = 0		#count matches
						refpsites = (int(st),int(st)+len(lp)-1)
					#	refp = ref[int(st)-1:int(st)+len(lp)-1]	#peptide in ref
						refp = vaccines[int(st)-1:int(st)+len(lp)-1] #peptide in ref
						for i in range(len(lp)):
							if lp[i]== refp[i]:
								m = m+1
							else:
								pass
						if float(m)/len(lp)>=0.7:
						#	print st, float(m)/len(lp), lp, refp			
							if hh not in matches:
								matches[hh] = [[ph],[st]]
							else:
								if st not in matches[hh][-1]:
									matches[hh][-1].append(st)
								if ph not in matches[hh][0]:
									matches[hh][0].append(ph)
							if ph not in Pmatches:
								Pmatches[ph] = range(int(st)-1,len(lp)+int(st)-1)
								if hh[0]=="A" or hh[0]=="B":
								#	print hh
									AB[ph] = range(int(st)-1,len(lp)+int(st)-1)
							else:
								pss = range(int(st)-1,len(lp)+int(st)-1)
								for ps in pss:
									if ps not in Pmatches[ph]:
										Pmatches[ph].append(ps)
								if hh[0]=="A" or hh[0]=="B":
									for ps in pss:
										if ps not in AB[ph]:
											AB[ph].append(ps)	
							
			match70(lanlhlas, all_matches, Pall_matches, Pall_matchesAB)
			match70(opt_lanlhlas, opt_matches, Popt_matches, Popt_matchesAB)

	#print Pall_matchesAB
	#get distances of autologous sequences from predicted reference sequence sites
        pidOUTall = {}                             #PID:mean dist all
	pidOUTall_AB  = {}
        pidOUTopt = {}                			#for optomals
	pidOUTopt_AB = {}
#	print Pall_matches
	def getdist(predicted, predOUT, preAB, preABout):			#Pall_matches,pidOUTall, Pall_matchesAB, pidOUTall_AB
		for pp in predicted:
			if pp in data:			#check if pid has sequence
				dist = 0
				for ppd in predicted[pp]:
					ep = int(ppd)
			
					dist = dist+int(matrix[vaccines[ep]][aacodes[data[pp][ep]]])
#					print vaccines[ep], data[pp][ep]
	#		print dist, len(predicted[pp])
				mdist = float("%0.1f"%(float(dist)/len(predicted[pp]))) #pid mean
				predOUT[pp] = str(mdist)
#				print pp, mdist
		#DO FOR OPTIMALS
		for po in preAB:
			if po in data:
				disto = 0
				for ppo in preAB[po]:
					epo = int(ppo)
					disto = disto+int(matrix[vaccines[epo]][aacodes[data[po][epo]]])
				mdisto = float("%0.1f"%(float(disto)/len(preAB[po]))) #pid mean
				preABout[po] = str(mdisto)
	getdist(Pall_matches,pidOUTall,Pall_matchesAB, pidOUTall_AB)
	getdist(Popt_matches,pidOUTopt, Popt_matchesAB, pidOUTopt_AB)	
	
	return (pidOUTall, pidOUTopt, pidOUTall_AB, pidOUTopt_AB, all_matches, opt_matches )			

predata = matchLANL(seqfile, hlafile, gene)
#print predata[0]
#outf1 = open(seqfile[:-4]+"_meanDistanceOUT.txt","w")
outf1 = open(outfile1,"w")
outf1.write("NB* The HLA-epitope data used for prediction is from 'http://www.hiv.lanl.gov/content/immunology/index.html' 2011 update\n\n-------------------------------------------\n\n")
outf1.write("Mean autologous peptide distances at ALL REACTIVE predicted epitope sites\n\nPID\tMeanSimilarity\n")
for pd in predata[0]:
	outf1.write(pd+"\t"+predata[0][pd]+"\n")
outf1.write("\n\n")
outf1.write("Mean autologous peptide distances at OPTIMAL predicted epitope sites\n\nPID\tMeanSimilarity\n")
for pd in predata[1]:
        outf1.write(pd+"\t"+predata[1][pd]+"\n")
outf1.write("\n\t------------------------------------\nRESULTS FOR HLA A and B alleles alone\n\nMean autologous peptide distances at ALL REACTIVE predicted HLA A & B epitope sites\n\nPID\tMeanSimilarity\n")
for aab in predata[2]:
	outf1.write(aab+"\t"+predata[2][aab]+"\n")
outf1.write("\n\nMean autologous peptide distances at OPTIMAL predicted HLA A & B epitope sites\n\nPID\tMeanSimilarity\n")
for oab in predata[3]:
        outf1.write(oab+"\t"+predata[3][oab]+"\n")
outf1.close()

#outf2 = open(seqfile[:-4]+"_PredictedSites.txt","w")
outf2 = open(outfile2,"w")
outf2.write("NB* The HLA-epitope data used for prediction is from 'http://www.hiv.lanl.gov/content/immunology/index.html' 2011 update\n\n\t-------------------------------------------\n\n")
outf2.write("All reactive Predicted epitope sites for each HLA\n\nHLA\tALLpredictions\tPIDSforALLpredictions\n")
for h in  predata[-2]:
	if len(predata[-2][h])>0:
		outf2.write(h+"\t"+str(predata[-2][h][-1])+"\t"+str(predata[-2][h][0])+"\n")

outf2.write("\n\n\t--------------------------------------------\n\nOPTIMAL Predicted epitope sites for each HLA\n\nHLA\tOptimalpredictions\tPIDSforOptimals\n")
for oh in  predata[-1]:
	if len(predata[-1][h])>0:
		outf2.write(oh+"\t"+str(predata[-1][oh][-1])+"\t"+str(predata[-1][oh][0])+"\n")
	
outf2.close()

