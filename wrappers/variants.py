#!/usr/bin/python
#get and list number of variants per amino acid position where differences exist between sequence and a reference sequence
#input is a tab_delimited alignment with the reference as the first sequence
#output is a text file listing each reference codon that has at least one variant and the frequencies of the variants
#ie. codon wt_residue No_variants Variant_details (AA-count;AA-count; etc)
#to run: python /home/nobubelo/uct_cbio_tools/src/tools/hiv_diversity_group/variants.py /home/nobubelo/uct_cbio_tools/src/tools/hiv_diversity_group/variant_alignment.txt /home/nobubelo/scratch/tests/varianttest.txt

import sys

infile, outfile = sys.argv[-2], sys.argv[-1]

file = open(infile,"r")
lines = file.readlines()

refid = lines[0].strip().split("\t")[0].strip()        #ref seq data
refseq = lines[0].strip().split("\t")[-1].strip()

variant = {}        #codon: [wt_residue, (Var1,countVar1), (Var2, countVar2)]

seqs = []        #all other sequences
for line in lines[1:]:                        #rest of the seqs
    id = line.strip().split("\t")[0].strip()
    seq = line.strip().split("\t")[-1].strip()
    if seq not in seqs:
        seqs.append(seq)

#do comparisons to refseq

for rf in range(len(refseq)):        #for each position
    wtype = refseq[rf]
    vars = {}            #list of variants for this position
    for sq in seqs:
        if sq[rf]!= wtype:
            if sq[rf] not in vars:
                vars[sq[rf]] = 1
            else:
                vars[sq[rf]] = vars[sq[rf]]+1
        else:
            pass
    #consolidate all into variant dict
    if len(vars)>0:
        variant[str(rf+1)] = [wtype, len(vars)]
        for v in vars:
            variant[str(rf+1)].append(v+"-"+str(vars[v]))

file.close()
outf = open(outfile,"w")
outf.write("codon\trefResidue\tVariantsCount\t; VariantsFrequency\n")
for vv in range(len(refseq)):
    v = str(vv+1)
    if v in variant:
        outf.write(v+'\t'+variant[v][0]+'\t'+str(variant[v][1])+'\t')
        for vr in variant[v][2:]:
            outf.write("; "+vr)
        outf.write("\n")
outf.close()
            