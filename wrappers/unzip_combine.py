#!/usr/bin/python
#uzip an archive and combine all file contents into one file
#to run eg: python unzip_combine.py jinnynew.zip jinnynew_ALLsequences.fas
import sys,os

infile, outf =  sys.argv[-2], sys.argv[-1]
cmd = "unzip "+infile +" -d "+infile[:-4]+"_filelist" 
os.system(cmd)

ndir = infile[:-4]+"_filelist"
st = os.listdir(ndir)
fpath = ndir+"/"
for sst in st:
	if os.path.isdir(fpath+sst):
		flist = os.walk(fpath+sst)
		for fl in flist:
			for i in range(len(fl[-1])):
				fpath = fl[0]+"/"+fl[-1][i]
				cmd = "cat "+fpath+" >> "+outf
				os.system(cmd)
	else:
		ffpath = ndir+"/"+sst
		cmd = "cat "+ffpath+" >> "+outf
		os.system(cmd)

rmc = "rm -R "+infile[:-4]+"_filelist"
os.system(rmc)
