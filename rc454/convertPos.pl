#! perl -w

# Copyright © 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE 
# This software and its documentation are the copyright of the 
# Broad Institute, Inc. All rights are reserved.
#  
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.
#
# This code was developed by Patrick Charlebois <patrickc@broadinstitute.org>

use strict;
use Getopt::Long;

my %option = (
	full	=> '',
	exact	=> '',
);

GetOptions(
  "full"	=> \$option{full},
  "exact"	=> \$option{exact},
) || die("Problem processing command-line options: $!\n");

my $poslist = shift || die;
my $reffile = shift || die;
my $inputfile = shift || die;
my $output = shift || die;

open(REFFILE, $reffile) || die;;
close REFFILE;
open(INPUTFILE, $inputfile) || die;
close INPUTFILE;
open(POSLIST, "$poslist") || die;;

system("cat $inputfile $reffile > $inputfile".".convmfa");

my $musclepath = "";
if($musclepath && substr($musclepath, length($musclepath) - 1) eq '/'){chop $musclepath;}
if($musclepath)
{
	system($musclepath."/muscle -in $inputfile".".convmfa -out $inputfile".".convafa");
}else{
	system("muscle -in $inputfile".".convmfa -out $inputfile".".convafa");
}

open(REFFILE, $reffile) || die;;
my $reffastaname = "";

while(my $line = <REFFILE>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		$reffastaname = $1;
		last;
	}
}

my $flag = 0;
my @features;

while(my $line = <POSLIST>)
{
	chomp($line);
	if($line =~ /\t/)
	{
		my @linearr = split(/\t/, $line);
		$features[$flag]{'name'} = $linearr[0];
		$features[$flag]{'start'} = $linearr[1];
		$features[$flag]{'end'} = $linearr[2];
		$flag++;
	}
}

open(AFA, "$inputfile".".convafa");

my $refseq = "";
my $otherseq = "";
my $curseq = "";

while (my $line = <AFA>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		if($curseq)
		{
			if($1 eq $reffastaname)
			{
				$otherseq = $curseq;
			}else{
				$refseq = $curseq;
			}
			$curseq = "";
		}
	}else
	{
		$curseq .= $line;
	}
}

if($refseq)
{
	$otherseq = $curseq;
}else{
	$refseq = $curseq;
}

my %posconversion;

my $refpos = 0;
my $otherpos = 0;

for(my $seqpos = 0; $seqpos < length($refseq); $seqpos++)
{
	my $curref = substr($refseq, $seqpos, 1);
	my $curother = substr($otherseq, $seqpos, 1);
	
	if($curref eq "-")
	{
		$otherpos++;
	}else{
		if($curother ne "-")
		{
			$otherpos++;
			$refpos++;
			$posconversion{$refpos} = $otherpos;
		}else{
			$refpos++;
			$posconversion{$refpos} = $otherpos."+";
		}
	}
}

my $lastotherpos = $otherpos;

open(OUTPUT, ">$output");

foreach my $feat (@features)
{

	my $start = $posconversion{$$feat{'start'}};
	my $end = $posconversion{$$feat{'end'}};

	my $full = 1;
	my $exact = 1;

	if($start == 0)
	{
		$start = "BeforeStart";
		$full = 0;
	}elsif($start eq $lastotherpos."+")
	{
		$start = "AfterEnd";
		$full = 0;
	}
	if($start =~ /(.+)\+/){
#		print "Start $start\n";
#		$start = $1;
		$exact = 0;
	}
	
	if($end == 0)
	{
		$end = "BeforeStart";
		$full = 0;
	}elsif($end eq $lastotherpos."+")
	{
		$end = "AfterEnd";
		$full = 0;
	}
	if($end =~ /(.+)\+/){
#		print "End $end\n";
#		$end = $1;
		$exact = 0;
	}
	
	
	if(($option{full} && $full == 0) || ($option{exact} && $exact == 0))
	{
		next;
	}
	print OUTPUT $$feat{'name'}."\t";
	print OUTPUT $start."\t".$end."\n";
}

system("rm $inputfile".".convmfa");
#system("rm $inputfile".".convafa");
