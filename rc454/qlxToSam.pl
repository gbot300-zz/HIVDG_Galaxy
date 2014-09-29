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
	bam	=> '',
	readfa => '',
	readq => '',
);

GetOptions(
  "bam"	=> \$option{bam},
  "readfa=s"	=> \$option{readfa},
  "readq=s"	=> \$option{readq},
);

my $qlxinput = shift || die("Usage : perl qlxToSam.pl input.qlx reference.fa output.sam [-bam]\n-bam will create a sorted bam output, but requires Newbler (use Newbler)\n");
my $ref = shift || die("Usage : perl qlxToSam.pl input.qlx reference.fa output.sam [-bam]\n-bam will create a sorted bam output, but requires Newbler (use Newbler)\n");
my $output = shift || die("Usage : perl qlxToSam.pl input.qlx reference.fa output.sam [-bam]\n-bam will create a sorted bam output, but requires Newbler (use Newbler)\n");

open(QLXFILE, $qlxinput) || die("Unable to open $qlxinput\n");
open(REF, $ref) || die("Unable to open $ref\n");
open(OUTPUT, ">$output") || die;

my $samtoolspath = "";
if($samtoolspath && substr($samtoolspath, length($samtoolspath) - 1) ne '/'){$samtoolspath .= '/';}

my %refseqs;
my $currefname = '';
my $currefseq = '';
while(my $line = <REF>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		if($currefseq)
		{
			$refseqs{$currefname} = $currefseq;
		}
		$currefname = $1;
		next;
	}else{
		$currefseq .= $line;
	}
}
$refseqs{$currefname} = $currefseq;

print OUTPUT "\@HD\tVN:1.0\tSO:coordinate\n";
foreach my $ref (sort keys %refseqs)
{
	print OUTPUT "\@SQ\tSN:".$ref."\tLN:".length($refseqs{$ref})."\n";
}

my $readid;
my $readstart;
my $readend;
my $readlen;
my $strand;
my $refstart;
my $refend;

my $refstr = '';
my $readstr = '';
my $qualstr = '';
my $mode;

my %reads;
my %quals;
my $curid = '';
if($option{readfa} && $option{readq})
{
	open(READS, $option{readfa});
	while(my $line = <READS>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$curid = $1;
		}else{
			$reads{$curid} .= $line;
		}
	}
	
	open(QUALS, $option{readq});
	while(my $line = <QUALS>)
	{
		chomp $line;
		if($line =~ />(.+)/)
		{
			$curid = $1;
		}else{
			if($line)
			{
				$quals{$curid} .= convertQual($line);
			}
		}
	}
}

while (my $line = <QLXFILE>)
{
	chomp $line;
#	print $line;
	if($line =~ />Read/)
	{
		my @linedata = split(/ /,$line);
		$readid = $linedata[1];
		$readstart = $linedata[2];
		$readend = $linedata[3];
		$readlen = $linedata[4];
		$strand = $linedata[5];
		$ref = $linedata[6];
		$refstart = $linedata[7];
		$refend = $linedata[8];
		$mode = 1;
	}elsif($mode == 1)
	{
		$refstr = $line;
		$mode = 2;
	}elsif($mode == 2)
	{
		$readstr = $line;
		$mode = 3;
	}elsif($mode == 3)
	{
		$mode = 4;
	}elsif($mode == 4)
	{
		$qualstr = $line;
		printSam();
		$mode = 0;
	}
}

if($option{bam})
{
	system($samtoolspath."samtools view -bS $output > $output.bam");
	system($samtoolspath."samtools sort $output.bam $output.sorted");
	system($samtoolspath."samtools index $output.sorted.bam");
}

sub printSam
{
	(my $cigarstr, my $parsedreadstr) = buildCigarStr($readstr, $refstr, $readstart, $readend, $readlen);
	print OUTPUT $readid."\t";
	if($strand eq '+')
	{
		print OUTPUT "0\t";
	}else{
		print OUTPUT "16\t";
		if($option{readfa} && $option{readq})
		{
			$reads{$readid} = revComp($reads{$readid});
			$quals{$readid} = reverse $quals{$readid};
		}
	}
	
	$qualstr =~ s/ //g;
	if($option{readfa} && $option{readq})
	{
		print OUTPUT "$ref\t$refstart\t255\t$cigarstr\t*\t0\t0\t".$reads{$readid}."\t".$quals{$readid}."\n";
	}else{
		print OUTPUT "$ref\t$refstart\t255\t$cigarstr\t*\t0\t0\t$parsedreadstr\t$qualstr\n";
	}
}

sub buildCigarStr
{
	my $curreadstr = shift;
	my $currefstr = shift;
	my $readstart = shift;
	my $readstop = shift;
	my $readlen = shift;
	
	my $cigarstr = '';
	
	my $count = 0;
	my $state = '';
	if($readstart > 1){
		$cigarstr .= ($readstart - 1)."S";
	}
	for(my $pos = 0; $pos < length($curreadstr); $pos++)
	{
		my $curread = substr($curreadstr, $pos, 1);
		my $curref = substr($currefstr, $pos, 1);
		if($curread eq '-')
		{
			if($state eq 'D'){$count++;}
			else{
				unless($state eq '')
				{
					$cigarstr .= $count.$state;
				}
				$count = 1;
				$state = 'D';
			}
		}elsif($curref eq '-')
		{
			if($state eq 'I'){$count++;}
			else{
				unless($state eq '')
				{
					$cigarstr .= $count.$state;
				}
				$count = 1;
				$state = 'I';
			}
		}else
		{
			if($state eq 'M'){$count++;}
			else{
				unless($state eq '')
				{
					$cigarstr .= $count.$state;
				}
				$count = 1;
				$state = 'M';
			}
		}
	}
	$cigarstr .= $count.$state;
	if($readstop < $readlen){
		$cigarstr .= ($readlen - $readstop)."S";
	}
	$curreadstr =~ s/\-//g;
	return $cigarstr, $curreadstr;
}

sub convertQual
{
	my $qualstr = shift;
	my @qualarr = split(/ /, $qualstr);
	my $ascQ = '';
	while(@qualarr)
	{
		$ascQ .= chr(33 + shift(@qualarr));
	}
	return $ascQ;
}

sub revComp
{
	my $seq = shift;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	my $rev = reverse $seq;
	return $rev;
}
