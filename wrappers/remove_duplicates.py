#!/usr/bin/python

#to remove duplicate sequences and renames identical seqnames:
#reads a tab_del alignment and outputs a fasta file with one copy of each sequence
#and prints deleted copies of duplicate seqs 
#also prints number of sequences in the output file

#TO RUN: python remove_duplicates.py duplicate_alignment.txt output.fasta log.txt

import sys

infile = sys.argv[-3]            #aafile, nucfile
outfile = sys.argv[-2]
logfile = sys.argv[-1]

singlecopies = []                       #list of non-duplicated sequences
data = {}                               #one_copy name:sequence
duplicates = []                         #list of duplicate copies
file = open(infile,"r")
lines = file.readlines()
for l in range(len(lines)):
    line = lines[l]
    name = line.strip().split("\t")[0]
    seq = line.strip().split("\t")[-1]
    if seq not in singlecopies:
        singlecopies.append(seq)
        if name not in data:    #avoid identical names for diff seqs
            data[name]=seq
        else:
#            print "renamed identical name for different sequences i.e., ", name
            data[name+"_r"+str(l)]= seq    
    else:
        duplicates.append(name)

outf = open(outfile,"w")
dfile = open(logfile,"w")
dfile.write("Duplicate sequences removed are:\n")

for d in data:
    outf.write(">"+d+"\n"+data[d]+"\n")

for dup in duplicates:
#    print "duplicate sequence removed is ", dup
    dfile.write(dup+"\n")

#print "number of remaining sequences in edited alignment is ",len(data)
dfile.write("\nNumber of remaining sequences in edited alignment is "+str(len(data)))

dfile.close()
file.close()
outf.close()