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
# This code was developed by Patrick Charlebois <patrickc@broadinstitute.org>

use strict;
use Getopt::Long;

my %option = (
	nqsmainqual     => 20,
	nqsareaqual     => 15,
	nqssize   		=> 5,
	nqvalue			=> 1,
);


GetOptions(
  "nqsmainqual=i"	=> \$option{nqsmainqual},
  "nqsareaqual=i"	=> \$option{nqsareaqual},
  "nqssize=i"		=> \$option{nqssize},
  "nqvalue=i"		=> \$option{nqvalue},
) || die("Problem processing command-line options: $!\n");

my $saminput = shift;
my $read = shift;
my $readqual = shift;
my $refinput = shift;
my $output = shift;

open(samFILE, $saminput) || die;
open(REF, $refinput) || die;
open(READ, $read) || die;
open(QUAL, $readqual) || die;
open(OUTPUT, ">$output") || die;

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

my $curreadid = '';
my %readseqs;
while(my $line = <READ>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		$curreadid = $1;
	}else{
		$readseqs{$curreadid} .= $line;
	}
}

my %readquals;
my %nqsquals;
while(my $line = <QUAL>)
{
	chomp $line;
	if($line =~ />(.+)/)
	{
		$curreadid = $1;
	}else{
		my @linequals = split(" ", $line);
		foreach my $linequal (@linequals)
		{
			$readquals{$curreadid} .= chr($linequal + 33);
		}
	}
}

foreach my $readid (keys %readquals)
{
	my $qualstr = $readquals{$readid};
	my $nqsstr = '';
	my $next5 = '';
	my $prev5 = '';
	for(my $qualpos = 0; $qualpos < length($qualstr); $qualpos++)
	{
		my $curqual = ord(substr($qualstr, $qualpos, 1)) - 33;
		$next5 = '';
		
		my $qualpos2 = $qualpos;
		while(length($next5) < $option{nqssize})
		{
			$qualpos2++;
			if($qualpos2 >= length($qualstr)){last;}
			if((ord(substr($qualstr, $qualpos2, 1)) - 33) == $option{nqvalue}){next;}
			
			if((ord(substr($qualstr, $qualpos2, 1)) - 33) < $option{nqsareaqual})
			{
				$next5 .= "0";
			}else{
				$next5 .= "1";
			}
		}
		if($curqual < $option{nqsmainqual})
		{
			$nqsstr .= "0";
		}elsif($next5 =~ /0/ || $prev5 =~ /0/)
		{
			$nqsstr .= "0";
		}else{
			$nqsstr .= "1";
		}
		if($curqual == $option{nqvalue}){next;}
		if($curqual >= $option{nqsareaqual}){$prev5 .= "1";}else{$prev5 .= "0";}
		if(length($prev5) > $option{nqssize}){$prev5 = substr($prev5, 1);}
	}
	$nqsquals{$readid} = $nqsstr;
}

my $readid;
my $readstart;
my $readend;
my $readlen;
my $strand;
my $ref;
my $refstart;
my $refend;
my $cigarstr = '';

my $refstr = '';
my $readstr = '';
my $qualstr = '';
my $mode;
my $rawreadstr;
my $rawqualstr;

while (my $line = <samFILE>)
{
	chomp $line;
	if($line =~ /^\@.*/)
	{
		next;
	}
		
	if(length($line) > 1)
	{
		my @linedata = split(/\t/,$line);
		$readid = $linedata[0];
		if($linedata[1] == 16)
		{
			$strand = '-'
		}else{
			$strand = '+'
		}
		$ref = $linedata[2];
		$refstart = $linedata[3];
		$cigarstr = $linedata[5];
		$rawreadstr = $linedata[9];
		$rawqualstr = $linedata[10];
		
		$readstr = '';
		$refstr = '';
		$qualstr = '';
		$refend = $refstart - 1;
		$readlen = 0;
		$readstart = 1;

		my $refpos = $refstart - 1;
		my $readpos = 0;
		while(length($cigarstr) > 0)
		{
			$cigarstr =~ /([0-9]+?)([A-Z]+?)(.*)/;
			my $len = $1;
			my $letter = $2;
			$cigarstr = $3;
		
			if($letter eq 'M')
			{
				$readlen += $len;
				$refend += $len;
				$readstr .= substr($rawreadstr, $readpos, $len);
				$qualstr .= substr($rawqualstr, $readpos, $len);
				$refstr .= substr($refseqs{$ref}, $refpos, $len);
				$readpos += $len;
				$refpos += $len;
			}elsif($letter eq 'D')
			{
				$refend += $len;
				$readstr .= "-" x $len;
				$qualstr .= " " x $len;
				$refstr .= substr($refseqs{$ref}, $refpos, $len);
				$refpos += $len;
			}elsif($letter eq 'I')
			{
				$readlen += $len;
				$refstr .= "-" x $len;
				$readstr .= substr($rawreadstr, $readpos, $len);
				$qualstr .= substr($rawqualstr, $readpos, $len);
				$readpos += $len;
			}
		}

		my $nqsstr = $nqsquals{$readid};
		$readend = $readstart + $readlen - 1;
		
		if($readlen != length($readseqs{$readid}))
		{
			if($strand eq '+')
			{
				$readseqs{$readid} =~ /(.*?)$rawreadstr(.*)/;
				$readstart = 1 + length($1);
				$readend = $readstart + $readlen - 1;
			}else{
				$readseqs{$readid} =~ tr/ACGTacgt/TGCAtgca/;
				my $revseq = reverse $readseqs{$readid};
				$revseq =~ /(.*?)$rawreadstr(.*)/;
				$readstart = 1 + length($2);
				$readend = $readstart + $readlen - 1;
			}
			$readlen = length($readseqs{$readid});
		}
		
		$nqsstr = substr($nqsstr,$readstart - 1, ($readend - $readstart + 1));
		
		if($strand eq '-'){$nqsstr = reverse $nqsstr;}

		my $fullnqsstr = '';
		my $nqspos = 0;
		for(my $i = 0; $i < length($readstr); $i++)
		{
			my $curnt = substr($readstr, $i, 1);
			if($curnt ne '-'){
				$fullnqsstr .= substr($nqsstr, $nqspos, 1);
				$nqspos++;
			}else{
				if(substr($fullnqsstr, length($fullnqsstr) - 1) == '0' || substr($fullnqsstr, length($fullnqsstr) - 1) == '9')
				{
					$fullnqsstr .= '9';
				}else{
					$fullnqsstr .= '8';
				}
			}
		}
		print OUTPUT ">Read $readid $readstart $readend $readlen $strand $ref $refstart $refend\n";
		print OUTPUT $refstr."\n".$readstr."\n".$fullnqsstr."\n".$qualstr."\n\n";
	}
}

