#!/usr/bin/python

'''

Author: Nobu

search for presence of  LANL epitopes (>=70% identity, for specific HLAs) in first sequence of an Amino Acid tab_del file eg 1st sequence could be a vaccine_insert for example, then determines divergence (using HIV-Between distance matrix) of autologous sequences from the 1st sequence at the identified epitope sites for each HLA allele genotype in the infected host: list of individual autologous epitope distances for each predicted epitope: the results can be used to compare these distances between 2 sets of data eg between placebo and vaccinees, using a non-parametric equivalent of a t-test

#USER INPUT = "tab del amino acid file (1st sequence assumed to be the reference) & "tab del HLA data file, no title line" & "gene name" which should be either Gag, Nef or Pol
#e.g. TO RUN: python matchEpitopeLANL.py lanl_seqtest.txt lanl_pHLA_all.txt Gag lanl_out1 lanl_out2 

#OUTPUT:- prints a list of Peptides, each peptide with list of distances of each autologus variant from submitted sequences  

Edited by Gordon Sep 2014
1. TODO Added Fasta input functionality
2. Wrote the transform_coords function to enable epitope searching in sequences that did not contain the entire HIV genome


'''

import sys, re, os

seqfile, hlafile, gene, outfile1 = sys.argv[-4], sys.argv[-3],sys.argv[-2], sys.argv[-1]
rel_start_pos = rel_end_pos = 0
#seqfile, hlafile, gene, outfile1, rel_start_pos, rel_end_pos = sys.argv[-6], sys.argv[-5],sys.argv[-4], sys.argv[-3], int(sys.argv[-2]), int(sys.argv[-1]) 		

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
	Function added by Gordon to transform the coords in the LANLallepitopes.txt file so that you can search for epitopes in your
	sequence, even if you only have a part of the HIV genome
	'''
	new_lines = []
	for line in lines:
		excl = False
		## extract the epitope positions from each line
		start, end = int(line.split('\t')[2]), int((line.split('\t')[3]))
		## transform positions 
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
# 	print data
	lanfile = open("/opt/exp_soft/galaxy-tools/custom_tools/meta_files/LANLallepitopes.txt", "r")			#LANLfile is either 'all' or 'alist' epitope file
	lanlines = lanfile.readlines()
	'''
	Gordon added the transform_coords function below to adjust the positional coordinates in the LANLallepitopes.txt file to enable 
	epitope searching if you dont have the entire HIV genome
	'''
	lanlines = transform_coords(lanlines)
	#print lanlines
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
				
#search for epitopes maching >= 70% on the ref seq for all PID hla alleles
	all_matches = {}			#hla:[[PID list][peptide_startsite list]]
	Pall_matches = {}			#peptide:[(pid,pept),(pid2,pept2), (),,,]
	refmatches = {}				#peptidesite:ref_petide sequence
	
	for ph in allphla:			#for each individual
		for hh in allphla[ph]:
			def match70(lanllist, matches, Pmatches):	#lanlhlas, all_matches Pallmatches, Pall_matchesAB OR opt_lanlhlas, opt_matches, Popt_matches, Popt_matchesAB
				if hh in lanllist:
		#			print hh
					for phd in lanllist[hh]:
						st = phd[0]		#start_site
						lp = phd[-1]		#LANL epitope sequence
						m = 0		#count matches
						refpsites = (int(st),int(st)+len(lp)-1)
					#	refp = ref[int(st)-1:int(st)+len(lp)-1]	#peptide in ref
						refp = vaccines[int(st)-1:int(st)+len(lp)-1] #peptide in ref
				#		print phd, st, lp
						for i in range(len(lp)):
							if lp[i]== refp[i]:
								m = m+1
							else:
								pass
						if float(m)/len(lp)>=0.7:
						#	print st, float(m)/len(lp), lp, refp			
							if hh not in matches:
								matches[hh] = [(ph,refpsites)]
							else:
								if (ph,refpsites) not in matches[hh]:
									matches[hh].append((ph,refpsites))
								else:
									pass
							if ph in data:
								pidpep = data[ph][int(st)-1:int(st)+len(lp)-1]
								if refpsites not in Pmatches:
									Pmatches[refpsites] = [(ph,pidpep)]
								else:
									if (ph,pidpep) not in Pmatches[refpsites]:
										Pmatches[refpsites].append((ph,pidpep))
							if refpsites not in refmatches:
								refmatches[refpsites] = refp
							else:
								pass
							
			match70(lanlhlas, all_matches, Pall_matches)
			
	#get distances of autologous sequences from predicted reference sequence sites
        pidOUTall = {}                             #Peptide:mean dist all patients with restricting HLA allele

	def getdist(predicted, predOUT):			#Pall_matches,pidOUTall
		for pps in predicted:
		#	print pps, predicted[pps]
			for pp in predicted[pps]:
				if pp[0] in data:			#check if pid has sequence
					dist = 0
					ppds = range(int(pps[0]),int(pps[-1])+1)
					for ppd in ppds:
						ep = int(ppd)
						dist = dist+int(matrix[vaccines[ep]][aacodes[data[pp[0]][ep]]])
				
					mdist = "%0.1f"%(float(dist)/len(ppds)) #peptide mean
					if pps not in predOUT:
						predOUT[pps] = [(pp[0],mdist)]
					else:
						predOUT[pps].append((pp[0],mdist))
				#print pps, cout, dist, float(mdist)/10
	
	getdist(Pall_matches,pidOUTall)
	return (pidOUTall, all_matches, refmatches,Pall_matches)			

predata = matchLANL(seqfile, hlafile, gene)


#print predata[0]
outf1 = open(outfile1,"w")
outf1.write("NB* The HLA-epitope data used for prediction is from 'http://www.hiv.lanl.gov/content/immunology/index.html' 2011 update\n\n-------------------------------------------\n\n")
outf1.write("Mean autologous peptide distances at ALL REACTIVE predicted epitope sites\n\n")
for pd in predata[0]:
#	print pd, predata[0][pd]
	if len(predata[0][pd])>0:
		outf1.write(str(pd[0])+predata[-2][pd]+str(pd[-1])+"\t")
		restrict=[]
		for hr in predata[1]:
			for hhr in predata[1][hr]:
				if pd in hhr:
					if hr not in restrict:
						restrict.append(hr)
					else:
						pass
		outf1.write(str(restrict)+"\n")

		for si in predata[0][pd]:
			outf1.write(si[0]+"\t")
			for pa in predata[-1][pd]:
				if pa[0]== si[0]:
					outf1.write(pa[-1]+"\t"+str(si[-1])+"\n")
		outf1.write("\n\n")
outf1.close()

