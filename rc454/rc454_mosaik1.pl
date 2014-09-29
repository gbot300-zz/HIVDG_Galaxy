#!/usr/bin/env perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE 
# This software and its documentation are the copyright of the 
# Broad Institute, Inc. All rights are reserved.
#  
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.
#
# This code was developed by Patrick Charlebois <patrickc@broadinstitute.org> with support from
# Matthew Henn <mhenn@broadinstitute.org>

use strict;
use Getopt::Long;

my %option = (
	minhomosize		=> 2,
	nqvalue 		=> 1,
	nocafie			=> 0,
	nqsmainqual     => 20,
	nqsareaqual     => 15,
	nqssize   		=> 5,
	details			=> 0,
	slicesize		=> 5,
	gap3window		=> 10,
	primers			=> '',
	minpctvar		=> 0.25,
	minnbvar		=> 5,
	genelist 		=> '',
	noorf 			=> '',
	bam				=> '',
	noclean			=> '',
	primerbuffer	=> 2,
);
GetOptions(
  "nocafie"				=> \$option{nocafie},
  "minhomosize=i"		=> \$option{minhomosize},
  "nqvalue=i"			=> \$option{nqvalue},
  "nqsmainqual=i"		=> \$option{nqsmainqual},
  "nqsareaqual=i"		=> \$option{nqsareaqual},
  "nqssize=i"			=> \$option{nqssize},
  "details"				=> \$option{details},
  "slicesize=i"			=> \$option{slicesize},
  "gap3window=i"		=> \$option{gap3window},
  "primers=s"			=> \$option{primers},
  "minpctvar=f"			=> \$option{minpctvar},
  "minnbvar=i"			=> \$option{minnbvar},
  "genelist=s"			=> \$option{genelist},
  "noorf"				=> \$option{noorf},
  "bam"					=> \$option{bam},
  "noclean"				=> \$option{noclean},
  "primerbuffer=i"		=> \$option{primerbuffer},
) || die("Problem processing command-line options: $!\n");



my $input = shift;
my $readfile = shift;
my $qualfile = shift;
my $reffile = shift;
my $output = shift;

my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
if($scriptpath && substr($scriptpath, length($scriptpath) - 1) eq '/'){chop $scriptpath;}

#Variables to read in .qln data
my $readid;
my $readstart;
my $readend;
my $readlen;
my $refstart;
my $refend;
my $strand;
my $contig;

my $minpos = 100000000000;
my $maxpos = 0;
my $nqsmainqual = 20;
my $mode = 0;

my %fixedie;
my $iecount = 0;

my %gappat;

#Hashes for problematic reads
my %unreliable;
my %probraw;
my %probfinal;

if($option{details})
{
	open(CHANGES, ">$output"."_changes.txt");
}
#### READ GENE LIST ####

my %primers;

if($option{primers})
{
	open(PRIMERS, $option{primers}) || die("Wrong Primer File\n");

	while(my $line = <PRIMERS>)
	{
		chomp($line);
		if($line =~ /\t/)
		{
			my @data = split(/\t/, $line);
			my $primori = '';
			if($data[0] =~ /.+\_F$/)
			{
				for(my $curprimpos = $data[1] - $option{primerbuffer}; $curprimpos <= $data[2]; $curprimpos++)
				{
					$primers{'F'}{$curprimpos} = 1;
				}
			}elsif($data[0] =~ /.+\_R$/)
			{
				for(my $curprimpos = $data[1]; $curprimpos <= $data[2] + $option{primerbuffer}; $curprimpos++)
				{
					$primers{'R'}{$curprimpos} = 1;
				}
			}else{
				die("Primer names must end by _F if on the forward strand or _R if on the reverse strand");
			}
		}
	}
}

my %genepositions;

if($option{genelist})
{
	open(GENELIST, $option{genelist}) || die("Invalid genelist file\n");
	while(my $line = <GENELIST>)
	{
		chomp $line;
		if($line =~ /(.+?)\t(.+?)\t(.+)/)
		{
			for(my $i = $2; $i <= $3; $i++)
			{
				$genepositions{$i} = 1;
			}
		}
	}
}


my $refstring;
my $readstring;
my $qualstring;

my $firstline;
$readid = "";

my %quals;

open(QUAL, $qualfile);

print "Reading Qual file...\n";
while (my $line = <QUAL>)
{
	chomp $line;
	if($line =~ />(.*)/)
	{
		$readid = $1;
		$firstline = 1;
	}else{
		if($firstline != 1)
		{
			$quals{$readid} .= " ";
		}
		$quals{$readid} .= $line;
		$firstline = 0;
	}
}

open(READFILE, $readfile);
my $curid = "";
my %readseqs;
my %homoreadseqs;
my %homoqualseqs;
my %iereadseqs;
my %iequalseqs;

my %alignslices;

print "Reading Read file...\n";
while(my $line = <READFILE>)
{
	chomp $line;
	if($line =~ />(.*)/)
	{
		$curid = $1;
	}elsif($curid)
	{
		$readseqs{$curid} .= $line;
	}
}
close READFILE;

my %alignedids;
my $tempqln = "$output"."_homo";
my $tempfasta = "$output"."_homo.fasta";
my $tempqual = "$output"."_homo.qual";

open(FASTAREJECTED, ">$output"."_rejected.fasta");
open(QUALREJECTED, ">$output"."_rejected.qual");
open(PROBFINAL, ">$output"."_finalproblems.txt");

my $curheader = "";

my $cleanread = "";
my $cleanref = "";
my $cleanqual = "";

my %aligndata;

for(my $turn = 1; $turn <= 3; $turn++)
{
	if($turn == 1 && $option{nocafie}){next;}
	if($turn == 1)
	{
		print "Cleaning Incomplete Extension data...\n";
		open(INPUT, "$input");
		open(FASTACLEAN, ">$output"."_ie.fasta");
		open(QUALCLEAN, ">$output"."_ie.qual");
	}elsif($turn == 2)
	{
		print "Cleaning Homopolymer data...\n";
		if($option{nocafie})
		{
			open(INPUT, "$input");
		}else{
			open(INPUT, "$output"."_ie.qlx");
		}
		open(FASTACLEAN, ">$tempfasta");
		open(QUALCLEAN, ">$tempqual");
	}else{
		print "Cleaning other remaining gap...\n";
		open(INPUT, "$tempqln.qlx");
		open(FASTACLEAN, ">$output"."_cleaned.fasta");
		open(QUALCLEAN, ">$output"."_cleaned.qual");
	}
	
	undef %aligndata;
	undef %alignslices;
	
	############ IMPORT DATA FROM THE QLX INPUT #############	
	
	while(my $line = <INPUT>)
	{
		chomp $line;
		if($line =~ />Read.*/)
		{
			$curheader = $line;
			my @data = split(/\s/, $line);
			$readid = $data[1];
			$readstart = $data[2];
			$readend = $data[3];
			$readlen = $data[4];
			$strand = $data[5];
			$contig = $data[6];
			$refstart = $data[7];
			$refend = $data[8];

#			print $readid."\n";

			$alignedids{$readid} = 1;
			
			$mode = 1;
		}elsif(length($line) > 1 && $mode == 1){
			$refstring = $line;
			$mode = 2;
		}elsif(length($line) > 1 && $mode == 2){
			$readstring = $line;
			$mode = 3;
		}elsif(length($line) > 1 && $mode == 3){
			$qualstring = $line;

			### TRIM PRIMERS
			
			if($option{primers})
			{
				my $trimlen = 0;
				my $checkrefpos = $refstart;
				my $readstartdiff = $readstart - 1;
				$checkrefpos -= $readstartdiff;
					
					
				while($primers{'F'}{$checkrefpos})
				{
					$checkrefpos++;
					$trimlen++;
				}
				$trimlen -= $readstartdiff;
					
				if($trimlen > 0)
				{
					my $reftrimmed = 0;
					my $aligntrimlen = 0;
					while($reftrimmed < $trimlen)
					{
						my $curref = substr($refstring, 0, 1);
						if($curref ne '-')
						{
							$reftrimmed++;
						}
						$refstring = substr($refstring, 1);
						$aligntrimlen++;
					}
					my $readtrimmed = 0;
					for(my $i = 0; $i < $aligntrimlen; $i++)
					{
						my $curread = substr($readstring, 0, 1);
						if($curread ne '-')
						{
							$readtrimmed++;
						}
						$readstring = substr($readstring, 1);
					}
					$qualstring = substr($qualstring, $aligntrimlen);
					$refstart += $trimlen;
					$readseqs{$readid} = substr($readseqs{$readid}, $readstartdiff);
					for(my $i = 1; $i <= $readstartdiff; $i++)
					{
						$quals{$readid} =~ /.+?\s(.+)/;
						$quals{$readid} = $1;
					}
					
					$readseqs{$readid} = substr($readseqs{$readid}, $readtrimmed);
					for(my $i = 1; $i <= $readtrimmed; $i++)
					{
						$quals{$readid} =~ /.+?\s(.+)/;
						$quals{$readid} = $1;
					}
					$readlen -= $readstartdiff;
					$readlen -= $readtrimmed;
					$readend -= $readstartdiff;
					$readend -= $readtrimmed;
					$readstart = 1;
				}

				$trimlen = 0;
				$checkrefpos = $refend;

				$readstartdiff = $readstart - 1;
				$checkrefpos += $readstartdiff;

				while($primers{'R'}{$checkrefpos})
				{
					$checkrefpos--;
					$trimlen++;
				}
				$trimlen -= $readstartdiff;
					
				if($trimlen > 0){
					my $reftrimmed = 0;
					my $aligntrimlen = 0;
					while($reftrimmed < $trimlen)
					{
						my $curref = substr($refstring, length($refstring) - 1, 1);
						if($curref ne '-')
						{
							$reftrimmed++;
						}
						$refstring = substr($refstring, 0, length($refstring) - 1);
						$aligntrimlen++;
					}
					my $readtrimmed = 0;
					for(my $i = 0; $i < $aligntrimlen; $i++)
					{
						my $curread = substr($readstring, length($readstring) - 1, 1);
						if($curread ne '-')
						{
							$readtrimmed++;
						}
						$readstring = substr($readstring, 0, length($readstring) - 1);
					}
					$qualstring = substr($qualstring, 0, length($qualstring)-$aligntrimlen);
					$refend -= $trimlen;
					
					$readseqs{$readid} = substr($readseqs{$readid}, $readstartdiff);
					for(my $i = 1; $i <= $readstartdiff; $i++)
					{
						$quals{$readid} =~ /.+?\s(.+)/;
						$quals{$readid} = $1;
					}

					$readseqs{$readid} = substr($readseqs{$readid}, $readtrimmed);
					for(my $i = 1; $i <= $readtrimmed; $i++)
					{
						$quals{$readid} =~ /.+?\s(.+)/;
						$quals{$readid} = $1;
					}
					$readlen -= $readstartdiff;
					$readlen -= $readtrimmed;
					$readend -= $readstartdiff;
					$readend -= $readtrimmed;
					$readstart = 1;
				}
			}

			$aligndata{$readid}{readstart} = $readstart;
			$aligndata{$readid}{readend} = $readend;
			$aligndata{$readid}{readlen} = $readlen;
			$aligndata{$readid}{strand} = $strand;
			$aligndata{$readid}{contig} = $contig;
			$aligndata{$readid}{refstart} = $refstart;
			$aligndata{$readid}{refend} = $refend;
			$aligndata{$readid}{readstring} = $readstring;
			$aligndata{$readid}{refstring} = $refstring;
			$aligndata{$readid}{qualstring} = $qualstring;
			$mode = 4;
		}
	}

	unless($turn == 2)
	{
		print "building slices size $option{slicesize}\n";
		readSlices($option{slicesize});
	}
	
	foreach my $curreadid(sort keys %aligndata)
	{
		$readid = $curreadid;
		
		my $curpos = $aligndata{$readid}{refstart};
		my $fullqualstr;

		if($turn == 1 || ($turn == 2 && $option{nocafie}))
		{
			$fullqualstr = $quals{$readid};	
		}elsif($turn==2){
			$fullqualstr = $iequalseqs{$readid};
		}else{
			$fullqualstr = $homoqualseqs{$readid};
		}
			
		#### IF IT DOESN'T START AT POS 0 OF THE READ, CUT THE QUALS
		if($aligndata{$readid}{readstart} != 1)
		{
			for (my $i = 0; $i < ($aligndata{$readid}{readstart} - 1); $i++)
			{
				$fullqualstr =~ s/(.*?)\s(.*)/$2/;
			}
		}
		
			
		if($aligndata{$readid}{readend} != $aligndata{$readid}{readlen})
		{
			for (my $i = $aligndata{$readid}{readend}; $i < $aligndata{$readid}{readlen}; $i++)
			{
				$fullqualstr =~ s/(.*)\s(.+)/$1/;
			}
		}

		#IF NEG STRAND REVERSE QUALS
		if($aligndata{$readid}{strand} eq "-")
		{
			$fullqualstr = revQual($fullqualstr);
		}

		if($turn == 1 || $turn == 2)
		{
			$cleanread = "";
			$cleanref = "";
			$cleanqual = "";
		}

		################# READ CLEANING ##############

		if($turn == 1)
		{
			### FIX INCOMPLETE EXTENSIONS
			if($aligndata{$readid}{strand} eq "-")
			{
				my $revreadstr = reverse $aligndata{$readid}{readstring};
				my $revrefstr = reverse $aligndata{$readid}{refstring};
				my $revqualstr = revQual($fullqualstr);
				($cleanread, $cleanref, $cleanqual) = cleanIncompleteExtension($revreadstr, $revrefstr, $revqualstr, $readid);
				$cleanread = reverse $cleanread;
				$cleanref = reverse $cleanref;
				$cleanqual = revQual($cleanqual);
			}else{
				($cleanread, $cleanref, $cleanqual) = cleanIncompleteExtension($aligndata{$readid}{readstring}, $aligndata{$readid}{refstring}, $fullqualstr, $readid);
			}
		}elsif($turn == 2)
		{
			### REMOVE HOMOPOLYMER MISCOUNTS
			($cleanread, $cleanref, $cleanqual) = cleanHomo($aligndata{$readid}{readstring}, $aligndata{$readid}{refstring}, $fullqualstr, $readid);
		}else{
			### REMOVE OTHER GAPS UNLESS THEY ARE IN AN AREA WHERE IT BECOMES A MULTIPLE OF 3 ###
			($cleanread, $cleanref, $cleanqual) = cleanGaps($aligndata{$readid}{readstring}, $aligndata{$readid}{refstring}, $fullqualstr, $readid);
		}

		############ OUTPUT PRINTING ################
	
		### REMOVE GAPS AND PUT SEQUENCES IN THE ORIGINAL READ ORDER TO PRINT READS AND QUALS FILE
		if($aligndata{$readid}{strand} eq "-")
		{
			$cleanread = revDna($cleanread);
			$cleanqual = revQual($cleanqual);
		}

		while($cleanread =~ /(.*)\-(.*)/)
		{
			$cleanread = $1.$2;
		}
		if($cleanqual =~ /^ (.+)/)
		{
			$cleanqual = $1;
		}
		if($cleanqual =~ /(.+) $/)
		{
			$cleanqual = $1;
		}

		### PRINT COMPLETE READS
		
		printComplete($cleanread, $cleanqual, $turn);
		
		$mode = 0;
	}
	
	close INPUT;
	close FASTACLEAN;
	close QUALCLEAN;

	my $opttoadd = '-nqvalue='.$option{nqvalue}.' -nqsmainqual='.$option{nqsmainqual}.' -nqsareaqual='.$option{nqsareaqual}.' -nqssize='.$option{nqssize};

	if($turn == 1){
		print "IE COUNT : $iecount\n";
		if($scriptpath)
		{
			system("perl $scriptpath/runMosaik.pl $output"."_ie.fasta $output"."_ie.qual $reffile $output"."_ie $opttoadd");
		}else{
			system("perl runMosaik.pl $output"."_ie.fasta $output"."_ie.qual $reffile $output"."_ie $opttoadd");
		}
	}elsif($turn == 2){
		if($scriptpath)
		{
			system("perl $scriptpath/runMosaik.pl $tempfasta $tempqual $reffile $tempqln -gop=30 $opttoadd");
		}else{
			system("perl runMosaik.pl $tempfasta $tempqual $reffile $tempqln -gop=30 $opttoadd");
		}
		if($option{noorf})
		{
			if($option{details})
			{
				printRejected();
			}
			cleanFiles();
		}
	}
}

if($option{details})
{
	printRejected();
}

my $opttoadd = '-nqvalue='.$option{nqvalue}.' -nqsmainqual='.$option{nqsmainqual}.' -nqsareaqual='.$option{nqsareaqual}.' -nqssize='.$option{nqssize};

if($scriptpath)
{
	system("perl $scriptpath/runMosaik.pl $output"."_cleaned.fasta $output"."_cleaned.qual $reffile $output"."_final -gop=30 $opttoadd");
}else{
	system("perl runMosaik.pl $output"."_cleaned.fasta $output"."_cleaned.qual $reffile $output"."_final -gop=30 $opttoadd");
}

cleanFiles();


###################################################################################
#################################### SUBS #########################################
###################################################################################


### REVERSE DNA SEQUENCE ###
sub revDna
{
	my $readstr = shift;
	
	my $torev = $readstr;
	$torev =~ tr/ACGTacgt/TGCAtgca/;
	$readstr = reverse $torev;	
	return $readstr;
}

### REVERSE QUAL STRING SEQUENCE ###
sub revQual
{
	my $qualstr = shift;

	my $torev = $qualstr;
	if($torev =~ /(.*)\s+$/){$torev = $1;}
	$qualstr = "";
	while($torev =~ s/(.*)\s(.*)/$1/)
	{
		$qualstr .= $2." ";
	}
	$qualstr .= $torev;
	return $qualstr;
}


### DETERMINE THE PATTERN OF GAPS IN AN ALIGNMENT (WHAT POSITION AND LENGTH) ###
sub gapPattern
{
	my $qlnfile = shift;
	
	my $qpmode = 0;
	my $qpreadid;
	my $qpreadstart;
	my $qpreadend;
	my $qpreadlen;
	my $qpstrand;
	my $qpcontig;
	my $qprefstart;
	my $qprefend;
	my $qpreadstring = "";
	my $qprefstring = "";
	my $qpqualstring = "";
	
	undef %gappat;
	
	open(GAPPATIN, $qlnfile);
	while(my $line = <GAPPATIN>)
	{
		chomp $line;

		if($line =~ />(.*)/)
		{
			my @data = split(/\s/, $line);
			$qpreadid = $data[1];
			$qpreadstart = $data[2];
			$qpreadend = $data[3];
			$qpreadlen = $data[4];
			$qpstrand = $data[5];
			$qpcontig = $data[6];
			$qprefstart = $data[7];
			$qprefend = $data[8];
			$qpmode = 1;
		}elsif(length($line) > 1 && $qpmode == 1){
			$qprefstring = $line;
			$qpmode = 2;
		}elsif(length($line) > 1 && $qpmode == 2){
			$qpreadstring = $line;
			$qpmode = 3;
		}elsif(length($line) > 1 && $qpmode == 3){
			$qpqualstring = $line;
			
			my $qppos = $qprefstart;
			
			while($qppos <= $qprefend)
			{
				my $curread = substr($qpreadstring, 0, 1);
				my $curref = substr($qprefstring, 0, 1);
				my $gaplen = 0;
				if($curref eq "-" || $curread eq "-")
				{
					$gaplen = 0;
					my $flag = 0;
					if($curref eq "-")
					{
						while(substr($qprefstring, $flag, 1) eq '-')
						{
							$flag++;
						}
						$gaplen += $flag;
						$gappat{ref}{$qppos - 1}{$gaplen}++;
					}else{
						while(substr($qpreadstring, $flag, 1) eq '-')
						{
							$flag++;
						}
						$gaplen += $flag;
						$gappat{read}{$qppos}{$gaplen}++;
						$qppos += $gaplen;
					}
					$qpreadstring = substr($qpreadstring, $gaplen);
					$qprefstring = substr($qprefstring, $gaplen);
				}else{
					$gappat{$qppos}{$gaplen}++;
					$qppos++;
					$qpreadstring = substr($qpreadstring, 1);
					$qprefstring = substr($qprefstring, 1);
				}
			}
			
			$qpmode = 0;
		}
	}
	close GAPPATIN;
	return;
}

sub readSlices
{
	my $slicelength = shift;

	foreach my $curreadid (sort keys %aligndata)
	{
		if($aligndata{$curreadid}{refend} - $aligndata{$curreadid}{refstart} < 2 * $slicelength)
		{
			next;
		}
		
		my $refstr = $aligndata{$curreadid}{refstring};
		my $readstr = $aligndata{$curreadid}{readstring};
		my $refpos = $aligndata{$curreadid}{refstart} + $slicelength;
		my $currefstring = "";
		my $curreadstring = "";
		
		my $currefpos = $aligndata{$curreadid}{refstart};
		while($currefpos <=  $refpos + $slicelength)
		{
			my $curref = substr($refstr, 0, 1);
			my $curread = substr($readstr, 0, 1);
			
			$currefstring .= $curref;
			$curreadstring .= $curread;
			
			if($curref ne "-"){$currefpos++;}
			
			$refstr = substr($refstr, 1);
			$readstr = substr($readstr, 1);
		}
		
		my $readslice = $curreadstring;
		$readslice =~ s/\-//g;
		
		$alignslices{'size'.$slicelength}{$refpos}{$readslice}++;
		$alignslices{'size'.$slicelength}{$refpos}{'coverage'}++;
		
		while($refpos < $aligndata{$curreadid}{refend} - $slicelength)
		{
			my $curref = substr($refstr, 0, 1);
			my $curread = substr($readstr, 0, 1);
			
			$currefstring .= $curref;
			$curreadstring .= $curread;
			
			if($curref ne "-"){
				$refpos++;
				my $initgaplen = 0;
				if($currefstring =~ /^[ATGCN](\-+).*/)
				{
					$initgaplen = length($1);
				}
				$currefstring = substr($currefstring, $initgaplen + 1);
				$curreadstring = substr($curreadstring, $initgaplen + 1);
			
				$readslice = $curreadstring;
				$readslice =~ s/\-//g;
				
				$alignslices{'size'.$slicelength}{$refpos}{$readslice}++;
				$alignslices{'size'.$slicelength}{$refpos}{'coverage'}++;
			}

			$refstr = substr($refstr, 1);
			$readstr = substr($readstr, 1);
		}
	}
}


### SEARCHES IF A READ SLICE IS FOUND ENOUGH IN THE GENOME TO CONSIDER THAT IT'S NOT AN ERROR ###
sub sliceLookup
{
	my $readid = shift;
	my $refpos = shift;
	my $slicelength = shift;
	my $comptype = shift;
	
	my $currefpos = $aligndata{$readid}{refstart};
	
	my $alignpos = 0;
	
	my $readslice = '';
	my $refslice = '';
	my $startedwriting = 0;
	while($currefpos <= $refpos + $slicelength)
	{
		
		if($currefpos >= $refpos - $slicelength)
		{
			if($startedwriting == 0 && substr($aligndata{$readid}{refstring}, $alignpos, 1) eq "-")
			{
				$alignpos++;
				next;
			}
			$readslice .= substr($aligndata{$readid}{readstring}, $alignpos, 1);
			$refslice .= substr($aligndata{$readid}{refstring}, $alignpos, 1);
			$startedwriting = 1;
		}
		
		if(substr($aligndata{$readid}{refstring}, $alignpos, 1) ne "-")
		{
			$currefpos++;
		}
		$alignpos++;
	}

	$readslice =~ s/\-//g;
	
	unless($alignslices{'size'.$slicelength}{$refpos}{$readslice})
	{
		print "Problem in $readid with read slice $readslice at position $refpos\n";
	}
	
	if($alignslices{'size'.$slicelength}{$refpos}{$readslice} >= minSliceLookup($refpos, $slicelength, $comptype))
	{
		print PROBFINAL "$readid was not corrected in step #"."$comptype at pos $refpos : 'error' found ".$alignslices{'size'.$slicelength}{$refpos}{$readslice}." times with coverage of ".$alignslices{'size'.$slicelength}{$refpos}{'coverage'}."\n";
		return 1;
	}else{
		return 0;
	}
}


sub minSliceLookup
{
	my $refpos = shift;
	my $slicelength = shift;
	my $comptype = shift;
	
	# For Carry Forward
	if($comptype == 1 || $comptype == 2)
	{
		if($alignslices{'size'.$option{slicesize}}{$refpos}{'coverage'} > ($option{minnbvar}/$option{minpctvar}))
		{
			return int($alignslices{'size'.$option{slicesize}}{$refpos}{'coverage'} * $option{minpctvar});
		}else{
			return $option{minnbvar};
		}
	}else
	{
		return 2;
	}
}

sub printGapPat
{
	foreach my $refpos(sort keys %{$gappat{ref}})
	{
		foreach my $refgaplen(sort keys %{$gappat{ref}{$refpos}})
		{
			print "Ref\t$refpos\t$refgaplen\t".$gappat{ref}{$refpos}{$refgaplen}."\n";	
		}
	}
	print "\n\n\n";
	foreach my $readpos(sort keys %{$gappat{read}})
	{
		foreach my $readgaplen(sort keys %{$gappat{read}{$readpos}})
		{
			print "Read\t$readpos\t$readgaplen\t".$gappat{read}{$readpos}{$readgaplen}."\n";	
		}
	}
	return;
}

### REMOVE ALL HOMOPOLYMER GAPS UNLESS MULTIPLES OF 3 ###
sub cleanHomo
{
	my $readstr = shift;
	my $refstr = shift;
	my $qualstr = shift;
	my $curreadid = shift;

	my $clnread = "";
	my $clnref = "";
	my $clnqual = "";
	
	my $initqualstr = $qualstr;
	
	my $gaplen;
	my $curhomoreadstr;
	my $curhomoqualstr;
	my $currefpos = $aligndata{$curreadid}{refstart};
	my $curreadpos = $aligndata{$curreadid}{readstart};
	
	for(my $pos = 0; $pos < length($refstr); $pos++)
	{
		my $curread = substr($readstr, $pos, 1);
		my $curref = substr($refstr, $pos, 1);
		my $curqual;

		my $inhomo = 0;
		my $homobase = $curread;
		$inhomo = inHomo($pos, $readstr, $refstr);

		unless($inhomo)
		{
			$inhomo = inHomo($pos, $refstr, $readstr);
		}
		
		if($curread ne "-")
		{
			if($qualstr =~ /(.+?) (.+)/)
			{
				$curqual = $1;
				$qualstr = $2;
			}else{
				$curqual = $qualstr;
			}
		}else{
			$curqual = $option{nqvalue};
		}
				
		if($curread eq "-" && $inhomo)
		{
			$gaplen = gapLength($readstr, $pos);
			if($gaplen % 3 != 0 || $gaplen == 3)
			{
				if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tInsert N\tHomopolymer\n";}
				$clnread .= "N";
				$clnref .= $curref;
				$clnqual .= "$curqual ";
			}else{
				$clnread .= "-";
				$clnref .= $curref;
			}
		}elsif($curref eq "-" && $inhomo)
		{
			$gaplen = gapLength($refstr, $pos);
			if($gaplen % 3 == 0)
			{
				$clnread .= $curread;
				$clnref .= $curref;
				$clnqual .= "$curqual ";
			}else{
				my $tempcurread = $curread;
				my $flagpos = $pos - 1;
				$nqsmainqual = $curqual;
				$minpos = $pos;
				my $nbgapplus = 0;
				my $nbgapminus = 0;
				
				while($flagpos >= 0)
				{
					$tempcurread = substr($readstr, $flagpos, 1);
					if($tempcurread eq '-')
					{
						$nbgapminus++;
						$flagpos--;
						next;
					}elsif($tempcurread ne $curread)
					{
						last;
					}
					my $tempcurqual = getQualAtPos($initqualstr, $readstr, $flagpos);
					if($tempcurqual < $nqsmainqual)
					{
						$nqsmainqual = $tempcurqual;
						$minpos = $flagpos;
					}
					$flagpos--;
				}
				$flagpos = $pos + 1;
				while($flagpos < length($readstr))
				{
					$tempcurread = substr($readstr, $flagpos, 1);
					if($tempcurread eq '-')
					{
						$nbgapplus++;
						$flagpos++;
						next;
					}elsif($tempcurread ne $curread)
					{
						last;
					}
					my $tempcurqual = getQualAtPos($initqualstr, $readstr, $flagpos);
					if($tempcurqual < $nqsmainqual)
					{
						$nqsmainqual = $tempcurqual;
						$minpos = $flagpos;
					}
					$flagpos++;
				}
				if($minpos < $pos)
				{
					my @clnqualarr = split(/ /, $clnqual);
					for(my $postoswap = (scalar(@clnqualarr) - ($pos - $minpos) + $nbgapminus); $postoswap < (scalar(@clnqualarr) - 1); $postoswap++)
					{
						$clnqualarr[$postoswap] = $clnqualarr[$postoswap + 1];
					}
					$clnqualarr[(scalar(@clnqualarr) - 1)] = $curqual;
					$clnqual = join(' ', @clnqualarr);
					$clnqual .= " ";
				}elsif($minpos > $pos){
					my @qualarr = split(/ /, $qualstr);
					if($minpos - $pos == 1)
					{
						my @qualarrrest = @qualarr[1..scalar(@qualarr)-1];
						@qualarr = ($curqual);
						push(@qualarr, @qualarrrest);
					}else{
						my @firstarr = @qualarr[0..($minpos-$pos-2)];
						my @secondarr = @qualarr[($minpos - $pos)..(scalar(@qualarr)-1)];
						@qualarr = ($curqual);
						push(@qualarr, @firstarr);
						push(@qualarr, @secondarr);
					}
					$qualstr = join(' ', @qualarr);
				}
				if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tRemoved $curread\tHomopolymer\n";}
			}
		}else{
			$clnread .= $curread;
			$clnref .= $curref;
			if($curread ne "-"){$clnqual .= "$curqual "};
		}
		if($curref ne '-'){$currefpos++;}
		if($curread ne '-'){$curreadpos++;}
	}
	return $clnread, $clnref, $clnqual;
}


### REMOVE ALL REMAINING GAPS UNLESS MULTIPLES OF 3 ###
sub cleanGaps
{
	my $readstr = shift;
	my $refstr = shift; 
	my $qualstr = shift;
	my $curreadid = shift;
	
	my $clnread = "";
	my $clnref = "";
	my $clnqual = "";
	
	my $gaplen = 0;
	my $readareagap = 0;
	my $refareagap = 0;
	my $gapremoved = 0;
	my $whileloopturns = 0;
	
	my $lastreadseq = '';
	my $lastrefseq = '';
#	print CHANGES $curreadid."\tINIT READ POS : ".$aligndata{$curreadid}{readstart}."\n";
	
	while(1)
	{
		my $refpos = 0;
		my $currefpos = $aligndata{$curreadid}{refstart};
		my $curreadpos = $aligndata{$curreadid}{readstart};
		for(my $pos = 0; $pos < length($refstr); $pos++)
		{
			my $updatedreadstring = $clnread.(substr($readstr, $pos));
			my $updatedrefstring = $clnref.(substr($refstr, $pos));
			my $updatedpos = length($clnread);
			
			my $curread = substr($readstr, $pos, 1);
			my $curref = substr($refstr, $pos, 1);
			my $curqual;
			
			if($curread ne "-")
			{
				if($qualstr =~ /(.+?) (.+)/)
				{
					$curqual = $1;
					$qualstr = $2;
				}else{
					$curqual = $qualstr;
				}
			}else{
				$curqual = $option{nqvalue};
			}
		
			if($curread eq "-" || $curref eq "-")
			{
				if($option{genelist})
				{
					my $fullrefpos = $aligndata{$curreadid}{refstart} + $refpos;
					unless($genepositions{$fullrefpos})
					{
#						print $fullrefpos."\n";
						$clnread .= $curread;
						$clnref .= $curref;
						if($curread ne "-")
						{
							$clnqual .= "$curqual ";
						}else{
							$refpos++;
						}
						next;
					}
				}
				my $gapdiff = abs(gapArea($updatedreadstring, $updatedpos) - gapArea($updatedrefstring, $updatedpos));

				my $slicelookuppos = ($aligndata{$curreadid}{refstart} + $refpos);
				if($aligndata{$curreadid}{refend} - $slicelookuppos < $option{slicesize})
				{
					$slicelookuppos =  $aligndata{$curreadid}{refend} - $option{slicesize};
				}elsif($slicelookuppos - $aligndata{$curreadid}{refstart} < $option{slicesize})
				{
					$slicelookuppos = $aligndata{$curreadid}{refstart} + $option{slicesize};
				}

				if($gapdiff % 3 != 0 && sliceLookup($curreadid, $slicelookuppos, $option{slicesize}, 2) == 0)
				{
					if($curread eq "-")
					{
						$gaplen = gapLength($updatedreadstring, $updatedpos);
						if($gaplen % 3 == 0)
						{
							$clnread .= $curread;
							$clnref .= $curref;
 						}else{
							$clnread .= "N";
							$clnref .= $curref;
							$clnqual .= "$curqual ";
							$gapremoved++;
							if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tAdd N\tNon-Homopolymer\n";}
						}
						$refpos++;
					}else{
						$gaplen = gapLength($updatedrefstring, $updatedpos);
						if($gaplen % 3 == 0)
						{
							$clnread .= $curread;
							$clnref .= $curref;
							$clnqual .= "$curqual ";
						}else{
							$gapremoved++;
							if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tRemoved $curread\tNon-Homopolymer\n";}
						}
					}
				}else{
					$clnread .= $curread;
					$clnref .= $curref;
					if($curread ne "-")
					{
						$clnqual .= "$curqual ";
					}else{
						$refpos++;
					}
				}
			}else{
				$clnread .= $curread;
				$clnref .= $curref;
				$clnqual .= "$curqual ";
				$refpos++;
			}
			if($curref ne '-'){
#				print CHANGES "Ref : ".$curreadid."\t".$currefpos."\t".$curref."\n";
				$currefpos++;
			}
			if($curread ne '-'){
#				print CHANGES "Read : ".$curreadid."\t".$curreadpos."\t".$curread."\n";
				$curreadpos++;
			}
		}
		my $refgapcount = $clnref =~ s/(\-)/$1/g;
		my $readgapcount = $clnread =~ s/(\-)/$1/g;
		$whileloopturns++;
		if(abs($refgapcount - $readgapcount) % 3 == 0 || $whileloopturns > 10)
		{
			if($whileloopturns > 10)
			{
				print $readid." didn't finish cleaning gaps after 10 turns\n";
			}
			last;
		}
		if($lastrefseq eq $clnref && $lastreadseq eq $clnread){last;}
		$readstr = $clnread;
		$refstr = $clnref;
		$qualstr = $clnqual;
		$clnread = "";
		$clnref = "";
		$clnqual = "";
		$lastrefseq = $refstr;
		$lastreadseq = $readstr;
	}

	if($gapremoved >= 3)
	{
		$unreliable{$readid} = 1;
		print PROBFINAL "$readid : removed $gapremoved gaps in non-homopolymer region, read unreliable\n";
	}

	return $clnread, $clnref, $clnqual;
}

sub cleanIncompleteExtension
{
	my $readstr = shift;
	my $refstr = shift;
	my $qualstr = shift;
	my $curreadid = shift;

	my $homobase = '';
	my $homolen = 0;

	my $readbaseno = 0;
	my $refbaseno = 0;

	my $initreadstr = $readstr;
	my $initrefstr = $refstr;
	my $currefpos = $aligndata{$curreadid}{refstart} - 1;
	my $curreadpos = $aligndata{$curreadid}{readstart} - 1;
	my $temprefpos = 0;
	my $tempreadpos = 0;

	for(my $curpos = 0; $curpos < length($refstr); $curpos++)
	{
		my $curread = substr($readstr, $curpos, 1);
		my $curref = substr($refstr, $curpos, 1);
		my $curqual;
		
		if($curread ne "-"){
			$readbaseno++;
			$curreadpos++;
		}
		if($curref ne "-"){
			$refbaseno++;
			$currefpos++;
		}

		if($curread eq $curref)
		{
			if($curread eq $homobase)
			{
				$homolen++;
			}else{
				$homolen = 1;
				$homobase = $curread;
			}
			next;
		}elsif($homolen < 2 || $curread eq $homobase || $curread eq "N" || $curref ne $homobase)
		{
			$homobase = '';
			$homolen = 0;
			next;
		}
		
		#MISMATCH FOUND AFTER HOMOPOLYMER
		my %countedbases;
		undef %countedbases;
		
		my %floworder;
		my $curfloworder = 1;
		my $setfloworder = 0;
		
		unless($curread eq '-')
		{
			$setfloworder = 1;
			$countedbases{$curread} = 1;
			if($curread eq 'T')
			{
				$floworder{T} = 1;
				$floworder{A} = 2;
				$floworder{C} = 3;
				$floworder{G} = 4;
			}elsif($curread eq 'A')
			{
				$floworder{T} = 4;
				$floworder{A} = 1;
				$floworder{C} = 2;
				$floworder{G} = 3;
			}elsif($curread eq 'C')
			{
				$floworder{T} = 3;
				$floworder{A} = 4;
				$floworder{C} = 1;
				$floworder{G} = 2;
			}elsif($curread eq 'G')
			{
				$floworder{T} = 2;
				$floworder{A} = 3;
				$floworder{C} = 4;
				$floworder{G} = 1;
			}
		}
		
		my $countdiffbases = 0;
		my $tempcurpos = $curpos + 1;
		my $tempreadpos = $curreadpos;
		my $temprefpos = $currefpos;
		if($tempcurpos >= length($refstr))
		{
			last;
		}
		my $lastbase = "";
		
		while($countdiffbases < 4)
		{
			my $tempcurread = substr($readstr, $tempcurpos, 1);
			my $tempcurref = substr($refstr, $tempcurpos, 1);
			if($tempcurread ne '-'){$tempreadpos++;}
			if($tempcurref ne '-'){$temprefpos++;}
			
			unless($setfloworder)
			{
				if($tempcurread ne '-' && $tempcurread ne 'N')
				{
					if($tempcurread eq 'T')
					{
						$floworder{T} = 1;
						$floworder{A} = 2;
						$floworder{C} = 3;
						$floworder{G} = 4;
					}elsif($tempcurread eq 'A')
					{
						$floworder{T} = 4;
						$floworder{A} = 1;
						$floworder{C} = 2;
						$floworder{G} = 3;
					}elsif($tempcurread eq 'C')
					{
						$floworder{T} = 3;
						$floworder{A} = 4;
						$floworder{C} = 1;
						$floworder{G} = 2;
					}elsif($tempcurread eq 'G')
					{
						$floworder{T} = 2;
						$floworder{A} = 3;
						$floworder{C} = 4;
						$floworder{G} = 1;
					}
					$curfloworder = 1;
					$setfloworder = 1;
				}
			}
			
			if($countedbases{$tempcurread} && $tempcurread ne $lastbase)
			{
				last;
			}
			if($tempcurread ne '-' && $tempcurread ne 'N' && $floworder{$tempcurread} < $curfloworder && $tempcurread ne $lastbase){
				last;
			} 
			unless($tempcurread eq '-' || $tempcurread eq 'N'){	$curfloworder = $floworder{$tempcurread};}
			
			if($tempcurread ne $tempcurref && $tempcurread eq $homobase)
			{
				#FOUND AN IE SUSPECT

				my $refpostoslice;
				if($aligndata{$readid}{strand} eq "-")
				{
					$refpostoslice = $aligndata{$readid}{refend} - $refbaseno;
					if($refbaseno < $option{slicesize}){$refpostoslice -= ($option{slicesize} - $refbaseno);}
					if($refpostoslice - $option{slicesize} < $aligndata{$readid}{refstart}){$refpostoslice += $option{slicesize} - ($refpostoslice - $aligndata{$readid}{refstart});}
				}else{
					$refpostoslice = $refbaseno + $aligndata{$readid}{refstart};
					if($refbaseno < $option{slicesize}){$refpostoslice += ($option{slicesize} - $refbaseno);}
					if($refpostoslice + $option{slicesize} > $aligndata{$readid}{refend}){$refpostoslice -= $option{slicesize} - ($aligndata{$readid}{refend} - $refpostoslice);}
				}

				if(sliceLookup($readid, $refpostoslice, $option{slicesize}, 1) == 1){last;}
			
				#FIX QUAL STRING
				$curqual = getQualAtPos($qualstr, $readstr, $tempcurpos);
				my $firstpartqual = getQualstrInRange($qualstr, $readstr, 0, $curpos - 1);
				my $secondpartqual = getQualstrInRange($qualstr, $readstr, $curpos, $tempcurpos - 1);
				my $thirdpartqual = getQualstrInRange($qualstr, $readstr, $tempcurpos + 1, length($refstr) - 1);
				
				$qualstr = $firstpartqual."  ".$curqual." ".$secondpartqual." ".$thirdpartqual;

				$qualstr =~ s/\s\s/ /g;

				#FIX READ STRING
				if($curread eq "-")
				{
					if($tempcurref eq "-"){
						$readstr = substr($readstr, 0, $curpos).$tempcurread.substr($readstr, $curpos + 1, ($tempcurpos-$curpos-1)).substr($readstr, $tempcurpos+1);
						$refstr = substr($refstr, 0, $tempcurpos).substr($refstr, $tempcurpos+1);
						if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tMoved $tempcurread from $tempreadpos and removed gap in ref\tCAFIE\n";}
					}else{
						$readstr = substr($readstr, 0, $curpos).$tempcurread.substr($readstr, $curpos + 1, ($tempcurpos-$curpos-1))."-".substr($readstr, $tempcurpos+1);
						if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tMoved $tempcurread from $tempreadpos and added gap in read\tCAFIE\n";}
					}
				}else{
					$readstr = substr($readstr, 0, $curpos).$tempcurread.substr($readstr, $curpos, ($tempcurpos-$curpos)).substr($readstr, $tempcurpos+1);
					if($option{details}){print CHANGES "$curreadid\t$curreadpos\t$currefpos\tMoved $tempcurread from $tempreadpos\tCAFIE\n";}
				}
				$iecount++;
				last;
			}else{
				if($tempcurread eq $homobase)
				{
					last;
				}elsif($tempcurread ne $lastbase && $tempcurread ne "-")
				{
					$countdiffbases++;
				}
				if($tempcurread ne "-")
				{
					$countedbases{$tempcurread} = 1;
				}
			}
			$tempcurpos++;
			$lastbase = $tempcurread;
			if($tempcurpos >= length($refstr)){last;}
		}
	}
	return $readstr, $refstr, $qualstr;
}

sub getQualSize
{
	my $qualstr = shift;
	my @qualarr = split(/ /, $qualstr);
	return scalar(@qualarr);
}


sub getQualAtPos
{
	my $qualstr = shift;
	my $readstr = shift;
	my $qualpos = shift; # 0 base
	
	my @qualarr = split(/\s/, $qualstr);
	
	my $cutreadstr = substr($readstr, 0, $qualpos);
	my $nbreadgap = ($cutreadstr =~ tr/\-/\./);
	if($qualpos >= scalar(@qualarr) + $nbreadgap)
	{
		print "ERROR getQualAtPos, $readid $qualpos >= ".scalar(@qualarr)."\n";
		die;
	}
	
	return $qualarr[$qualpos - $nbreadgap];
}

sub getQualstrInRange
{
	my $qualstr = shift;
	my $readstr = shift;
	my $firstpos = shift;
	my $lastpos = shift;

	my $firstreadstr = substr($readstr, 0, $firstpos);

	my $secondreadstr = substr($readstr, $firstpos, ($lastpos - $firstpos + 1));
	
	my $nbfirstreadgap = ($firstreadstr =~ tr/\-/\./);
	my $nbsecondreadgap = ($secondreadstr =~ tr/\-/\./);
	
	my $nbreadgap = $nbfirstreadgap + $nbsecondreadgap;
	
	my @qualarr = split(/ /, $qualstr);

	if($firstpos < 0 || $lastpos >= (scalar(@qualarr) + $nbreadgap))
	{
		print "ERROR getQualstrInRange $readid $nbreadgap $firstpos $lastpos ".scalar(@qualarr)."\n";
	}
	my $strtoreturn = "";
	for(my $curpos = $firstpos - $nbfirstreadgap; $curpos <= $lastpos - $nbreadgap; $curpos++)
	{
		if($curpos == ($firstpos - $nbfirstreadgap))
		{
			$strtoreturn = $qualarr[$curpos];
		}else{
			$strtoreturn .= " ".$qualarr[$curpos];
		}
	}
	
	return $strtoreturn;
}


### PRINT REJECTED FILES
sub printRejected
{
	open(READFILE, $readfile);
	my $printmode = 0;
	while (my $line = <READFILE>)
	{
		if ($line =~ />(.+)/)
		{
			unless($alignedids{$1})
			{
				print FASTAREJECTED ">$1\n";
				$printmode = 1;	
			}else{
				$printmode = 0;
			}
		}
		elsif($line =~ /(.+)/)
		{
			if($printmode == 1)
			{
				print FASTAREJECTED $line;
			}
		}
	}

	open(QUALFILE, $qualfile);
	$printmode = 0;
	while (my $line = <QUALFILE>)
	{
		if ($line =~ />(.+)/)
		{
			unless($alignedids{$1})
			{
					print QUALREJECTED ">$1\n";
				$printmode = 1;	
			}else{
				$printmode = 0;
			}
		}
		elsif($line =~ /(.+)/)
		{
			if($printmode == 1)
			{
				print QUALREJECTED $line;
			}
		}
	}
}

sub printComplete
{
	my $clnread = shift;
	my $clnqual = shift;
	my $turn = shift;
	

	if($aligndata{$readid}{readstart} != 1)
	{
		my $toadd;
		if($turn == 1 || ($turn==2 && $option{nocafie}))
		{
			$toadd = substr($readseqs{$readid}, 0, ($aligndata{$readid}{readstart} - 1));
		}elsif($turn==2)
		{
			$toadd = substr($iereadseqs{$readid}, 0, ($aligndata{$readid}{readstart} - 1));
		}else{
			$toadd = substr($homoreadseqs{$readid}, 0, ($aligndata{$readid}{readstart} - 1));
		}
		$clnread = $toadd.$clnread;
		
		my $extraqual = "";
		
		for(my $i = 0; $i < ($aligndata{$readid}{readstart} - 1); $i++)
		{
			if($turn == 1 || ($turn==2 && $option{nocafie}))
			{
				$quals{$readid} =~ s/(.*?)\s(.*)/$2/;
				$extraqual .= $1." ";
			}elsif($turn == 2){
				$iequalseqs{$readid} =~ s/(.*?)\s(.*)/$2/;
				$extraqual .= $1." ";
			}else{
				$homoqualseqs{$readid} =~ s/(.*?)\s(.*)/$2/;
				$extraqual .= $1." ";
			}
		}
		
		$clnqual = $extraqual.$clnqual;
	}
	
	if($aligndata{$readid}{readend} != $aligndata{$readid}{readlen})
	{
		my $diff = $aligndata{$readid}{readlen} - $aligndata{$readid}{readend};
		my $toadd;
		if($turn == 1 || ($turn==2 && $option{nocafie}))
		{
			$toadd = substr($readseqs{$readid}, (length($readseqs{$readid}) - $diff));
		}elsif($turn==2){
			$toadd = substr($iereadseqs{$readid}, (length($iereadseqs{$readid}) - $diff));
		}else{
			$toadd = substr($homoreadseqs{$readid}, (length($homoreadseqs{$readid}) - $diff));
		}
		
		$clnread .= $toadd;
		
		my $extraqual = "";
		for (my $i = 0; $i < $diff; $i++)
		{
			if($turn == 1 || ($turn==2 && $option{nocafie}))
			{
				$quals{$readid} =~ s/(.+)\s(.+)/$1/;
				$extraqual = $2." ".$extraqual;
			}elsif($turn == 2)
			{
				$iequalseqs{$readid} =~ s/(.+)\s(.+)/$1/;
				$extraqual = $2." ".$extraqual;
			}else{
				$homoqualseqs{$readid} =~ s/(.+)\s(.+)/$1/;
				$extraqual = $2." ".$extraqual;
			}
		}
		chop $extraqual; #remove trailing space
		$clnqual .= " ".$extraqual;
	}
	print FASTACLEAN ">".$readid."\n";
	print FASTACLEAN $clnread."\n";
	print QUALCLEAN ">".$readid."\n";
	
	if($clnqual =~ /\s$/)
	{
		chop $clnqual;
	}
	if($clnqual =~ /\s\s/)
	{
		$clnqual =~ s/\s\s/ /g;
	}
	
	print QUALCLEAN $clnqual."\n";
	if($turn == 1)
	{
		$iereadseqs{$readid} = $clnread;
		$iequalseqs{$readid} = $clnqual;
	}elsif($turn == 2)
	{
		$homoreadseqs{$readid} = $clnread;
		$homoqualseqs{$readid} = $clnqual;
	}
}

### DETERMINE IF A POSITION IS IN AN HOMOPOLYMER IN A SEQUENCE
sub inHomo
{
	my $postocheck = shift;
	my $dnastring = shift;
	my $compstring = shift;
	my $homostrres = 1;
	my $homores = substr($dnastring, $postocheck, 1);

	my $previousmatch = 0;

	for(my $homopos = $postocheck - 1; $homopos >= $postocheck - $option{minhomosize}; $homopos--)
	{
		if($homopos < 0){last;}

		my $resnow = substr($dnastring, $homopos, 1);
		my $compres = substr($compstring, $homopos, 1);

		if($resnow eq $homores)
		{
			if($compres eq $resnow){$previousmatch = 1;}
			$homostrres++;
		}else{last;}

		if($homostrres == $option{minhomosize} && $previousmatch)
		{
			return $homores;
		}
	}
	
	for(my $homopos = $postocheck + 1; $homopos <= $postocheck + $option{minhomosize}; $homopos++)
	{
		if($homopos >= length($dnastring)){last;}
		
		my $resnow = substr($dnastring, $homopos, 1);
		my $compres = substr($compstring, $homopos, 1);

		if($resnow eq $homores)
		{
			if($compres eq $resnow){$previousmatch = 1;}
			$homostrres++;
		}else{last;}
		
		if($homostrres == $option{minhomosize} && $previousmatch)
		{
			return $homores;
		}
	}
	return 0;
}

### CALCULATE HOW MANY GAPS THAT ARE NOT MULT.3 ARE IN THE AREA AROUND A CENTRAL GAP
sub gapArea
{
	my $dnastring = shift;
	my $postocheck = shift;
	my $totalgap = 0;
	my $areasize = $option{gap3window};

	my $curpos = $postocheck;
	
	my $initialgaplen = 0;
	my $initialposlow;
	my $initialposhigh;
	
	while(my $resnow = substr($dnastring, $curpos, 1))
	{
		if($resnow eq "-")
		{
			$initialgaplen++;
		}else{
			$initialposlow = $curpos;
			last;
		}
		$curpos--;
		if($curpos < 0)
		{
			last;
		}
	}
	
	$curpos = $postocheck + 1;
	while(my $resnow = substr($dnastring, $curpos, 1))
	{
		if($resnow eq "-")
		{
			$initialgaplen++;
		}else{
			$initialposhigh = $curpos;
			last;
		}
		$curpos++;
		if($curpos > length($dnastring))
		{
			last;
		}
	}

	if($initialgaplen % 3 != 0)
	{
		$totalgap += $initialgaplen;
	}
	
	my $curgaplen = 0;

	$curpos = $initialposlow;
	
	while(1)
	{
		if(substr($dnastring, $curpos, 1) eq "-")
		{
			$curgaplen++;
		}else{
			if($curgaplen)
			{
				if($curgaplen % 3 != 0)
				{
					$totalgap += $curgaplen;
				}
				$curgaplen = 0;
			}
			if(abs($postocheck - $curpos) > $areasize)
			{
				last;
			}
		}
		$curpos--;
		if($curpos < 0)
		{
			last;
		}
	}

	$curpos = $initialposhigh;
	$curgaplen = 0;
	while(1)
	{
		if(substr($dnastring, $curpos, 1) eq "-")
		{
			$curgaplen++;
		}else{
			if($curgaplen)
			{
				if($curgaplen % 3 != 0)
				{
					$totalgap += $curgaplen;
				}
				$curgaplen = 0;
			}
			if(abs($postocheck - $curpos) > $areasize)
			{
				last;
			}
		}
		$curpos++;
		if($curpos > length($dnastring))
		{
			last;
		}	
	}

	return $totalgap;
}


### CALCULATE THE LENGTH OF A GAP
sub gapLength
{
	my $dnastring = shift;
	my $postocheck = shift;

	my $resnow = "-";
	my $gaplength = 1;

	my $flag = 1;
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck - $flag), 1);
		if($resnow eq "-")
		{
			$gaplength++;
			$flag++;
		}
	}

	$flag = 1;
	$resnow = "-";
	while($resnow eq "-")
	{
		$resnow = substr($dnastring, ($postocheck + $flag), 1);
		if($resnow eq "-")
		{
			$gaplength++;
			$flag++;
		}
	}
	return $gaplength;
}

sub cleanFiles
{
	if($option{noorf})
	{
		system("mv $output"."_homo.bam $output"."_final.bam");
		system("mv $output"."_homo.fasta $output"."_cleaned.fasta");
		system("mv $output"."_homo.qual $output"."_cleaned.qual");
		system("mv $output"."_homo.qlx $output"."_final.qlx");
		system("mv $output"."_homo.sam $output"."_final.sam");
	}

	unless($option{details})
	{
		system("rm $output"."_finalproblems.txt");
		unless($option{noorf})
		{
			system("rm $output"."_homo.bam");
			system("rm $output"."_homo.fasta");
			system("rm $output"."_homo.qual");
			system("rm $output"."_homo.qlx");
			system("rm $output"."_homo.sam");
		}

		system("rm $output"."_ie.bam");
		system("rm $output"."_ie.fasta");
		system("rm $output"."_ie.qual");
		system("rm $output"."_ie.qlx");
		system("rm $output"."_ie.sam");
		
		system("rm $output"."_rejected.fasta");
		system("rm $output"."_rejected.qual");
	}
	
	unless($option{bam})
	{
		system("rm $output"."_raw.bam");
		system("rm $output"."_raw.sam");
		system("rm $output"."_final.bam");
		system("rm $output"."_final.sam");
	}
	
	if($option{noclean})
	{
		system("rm $output"."_cleaned.fasta");
		system("rm $output"."_cleaned.qual");
	}
	
	print "Done\n";
	exit;
}