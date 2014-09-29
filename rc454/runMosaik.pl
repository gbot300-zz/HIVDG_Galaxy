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
# This wrapper is designed to run the software Mosaik which can be found at 
# http://bioinformatics.bc.edu/marthlab/Mosaik

use strict;
use Getopt::Long;
use warnings;

my %option = (
	hs	=> 10,
	act	=> 15,
	mmp	=> 0.25,
	minp	=> 0.25,
	mms	=> -9,
	hgop => 4,
	gop	=> 15,
	gep	=> 6.66,
	nqsmainqual     => 20,
	nqsareaqual      => 15,
	nqssize   	=> 5,
	bw			=> 0,
	nqvalue	=> 1,
	unique 	=> 1,
	st	=> '454',
	qlxonly => 0,
);
GetOptions(
  "hs=i"	=> \$option{hs},
  "act=i"	=> \$option{act},
  "mmp=f"	=> \$option{mmp},
  "minp=f"	=> \$option{minp},
  "mms=i"	=> \$option{mms},
  "hgop=i"	=> \$option{hgop},
  "gop=i"	=> \$option{gop},
  "gep=i"	=> \$option{gep},
  "bw=i"	=> \$option{bw},
  "nqsmainqual=i"	=> \$option{nqsmainqual},
  "nqsareaqual=i"	=> \$option{nqsareaqual},
  "nqssize=i"	=> \$option{nqssize},
  "nqvalue=i"	=> \$option{nqvalue},
  "unique=i"	=> \$option{unique},
  "st=s"	=> \$option{st},
  "qlxonly" => \$option{qlxonly},
) || die("Problem processing command-line options: $!\n");

my $inputreads = shift;
my $inputqual = shift;
my $inputref = shift;
my $output = shift;

my $mosaikpath = "";
if($mosaikpath && substr($mosaikpath, length($mosaikpath) - 1) eq '/'){chop $mosaikpath;}

my $refdat = $output.".refdat";
my $readdat = $output.".readdat";

my $unique = "unique";
my $uniquesort = "-u";
if($option{unique} == 0)
{
	$unique = "all";
	$uniquesort = "";
}

my $bw = '';
if($option{bw})
{
	$bw = " -bw ".$option{bw};
}


#Build Ref and Reads .dat with Mosaik
if($mosaikpath)
{
	system($mosaikpath."/MosaikBuild -fr $inputref -oa $refdat");
	system($mosaikpath."/MosaikBuild -fr $inputreads -fq $inputqual -out $readdat -tn 500 -st ".$option{st});
	system($mosaikpath."/MosaikAligner -in $readdat -out $readdat.vsRef -ia $refdat -hs ".$option{hs}." -act ".$option{act}." -mm 500 -mmp ".$option{mmp}." -mmal -minp ".$option{minp}." -gop ".$option{gop}." -hgop ".$option{hgop}." -gep ".$option{gep}."$bw -m $unique");
	system($mosaikpath."/MosaikSort -in $readdat.vsRef -out $readdat.sorted");
	system($mosaikpath."/MosaikText -in $readdat.sorted -axt $output.axt $uniquesort");
}else{
	system("MosaikBuild -fr $inputref -oa $refdat");
	system("MosaikBuild -fr $inputreads -fq $inputqual -out $readdat -tn 500 -st ".$option{st});
	system("MosaikAligner -in $readdat -out $readdat.vsRef -ia $refdat -hs ".$option{hs}." -act ".$option{act}." -mm 500 -mmp ".$option{mmp}." -mmal -minp ".$option{minp}." -gop ".$option{gop}." -hgop ".$option{hgop}." -gep ".$option{gep}."$bw -m $unique");
	system("MosaikSort -in $readdat.vsRef -out $readdat.sorted");
	system("MosaikText -in $readdat.sorted -axt $output.axt $uniquesort");
}

unless($option{qlxonly})
{
	if($mosaikpath)
	{
		system($mosaikpath."/MosaikText -in $readdat.sorted -sam $output.sam $uniquesort");
		system($mosaikpath."/MosaikText -in $readdat.sorted -bam $output.bam $uniquesort");
	}else{
		system("MosaikText -in $readdat.sorted -sam $output.sam $uniquesort");
		system("MosaikText -in $readdat.sorted -bam $output.bam $uniquesort");
	}
}
system("gzip -d -f $output.sam.gz");

print "\nCreating .qlx file...\n";
my $toadd = " -nqvalue=".$option{nqvalue};
if($option{nqsmainqual} != 20)
{
	$toadd .= " -nqsmainqual=".$option{nqsmainqual};
}
if($option{nqsareaqual} != 15)
{
	$toadd .= " -nqsareaqual=".$option{nqsareaqual};
}
if($option{nqssize} != 5)
{
	$toadd .= " -nqssize=".$option{nqssize};
}


# Global variables for axtToQlx
my %quals;
my %nqsquals;
my $readstart;
my $readend;
my $contig;
my $readstring = "";
my $refstring = "";
my $lastpos;
my %ntcount;
my %nthqcount;
my %refseq;
my %strand;
my @query;
my $state = 0;
my $refstart;
my %readlen;


open (QLXOUT, ">$output.qlx") || die;
axtToQlx("$output.axt", $inputqual);

print "\nRemoving temp files...\n";
system("rm $readdat");
system("rm $refdat");
system("rm $readdat.vsRef");
system("rm $readdat.sorted");
system("rm $output.axt");

sub axtToQlx
{
	my $axtfile = shift;
	my $qualfile = shift;
	
	open (AXT, $axtfile) || die;
	open (QUAL, $qualfile) || die;
	
	print "NQS Parameters :\n";
	print "Minimum base quality : ".$option{nqsmainqual}."\n";
	print "Minimum neighborhood base quality : ".$option{nqsareaqual};
	print " (every base in the neighborhood)\n";
	print "Neighborhood size on each side of the base : ".$option{nqssize}."\n";
	
	#***** FIX BY SWITCHING QUALS TO ARRAY, THEN CAN ADD QLX EASY *****
	### Read QUALS

	my $nqsmainqual = $option{nqsmainqual};
	my $nqsareaqual = $option{nqsareaqual};
	my $nqrange = $option{nqssize};

	my $readid = "";
	
	print "Read Quals:\n";
	while (my $line = <QUAL>)
	{
		chomp $line;
		if($line =~ />(.*)/)
		{
			if($quals{$readid})
			{
				$readlen{$readid} = scalar(@{$quals{$readid}});
			}
			$readid = $1;
		}else{
			
			my @linedata = split(/ /, $line);
			foreach my $curqualdata (@linedata)
			{
				push(@{$quals{$readid}}, $curqualdata);
			}
		}
	}
	$readlen{$readid} = scalar(@{$quals{$readid}});
	
	### Eval QUALS
	print "\nEval Quals:\n";
	foreach my $read (keys %quals)
	{
		my @seqquals = @{$quals{$read}};
		my $qualflag = 0;
		
		foreach my $qscore(@seqquals)
		{
			if($qscore eq ''){print "problem $read\n";}
			if($qscore < $nqsmainqual)
			{
				push(@{$nqsquals{$read}}, "0");
				$qualflag++;
				next;
			}
			my $ok = 1;
				
			my $flag = $qualflag - 1;
			my $basecounted = 0;
			while($basecounted < $nqrange)
			{
				if($flag < 0)
				{
					last;
				}elsif($seqquals[$flag] == $option{nqvalue} && !($option{countN})){
					$flag--;
				}else{
					if($seqquals[$flag] < $nqsareaqual)
					{
						$ok = 0;
						last;
					}
					$flag--;
					$basecounted++;
				}
			}
			if($ok==1)
			{
				$flag = $qualflag + 1;
				$basecounted = 0;
				while($basecounted < $nqrange)
				{
					if($flag > (scalar(@seqquals) - 1))
					{
						last;
					}elsif($seqquals[$flag] == $option{nqvalue} && !($option{countN})){
						$flag++;
					}else{
						if($seqquals[$flag] < $nqsareaqual)
						{
							$ok = 0;
							last;
						}
						$flag++;
						$basecounted++;
					}
				}
			}
			if($ok == 1)
			{
				push(@{$nqsquals{$read}}, "1");
			}else{
				push(@{$nqsquals{$read}}, "0");
			}
			$qualflag++;
		}
	}
	
	$readid = "";

	### Read AXT output and assign codons
	while (my $line = <AXT>)
	{
		chomp $line;
		
		if($line =~ /[0-9]/)
		{
			if($readid)
			{
				evalRead($readstring, $refstring, $readid);
	#			print "OutEval\n";
			}
			@query = split(/\s/, $line);
			$contig = $query[1];
			$refstart = $query[2];
			$lastpos = $query[3];
			$readid = $query[4];
			$readstart = $query[5];
			$readend = $query[6];
			$strand{$readid} = $query[7];
			$readstring = "";
			$refstring = "";
			$state = 1;
			
			#### IF IT DOESN'T START AT POS 0 OF THE READ, CUT THE QUALS
			if($readstart > 1)
			{
				for(my $i = $readstart; $i > 1; $i--)
				{
					shift(@{$quals{$readid}});
					shift(@{$nqsquals{$readid}});
				}
				my $lenbefore = scalar(@{$quals{$readid}});
			}
	
	
			if($readend < $readlen{$readid} && $strand{$readid} eq '-')
			{
				for(my $i = $readend; $i < $readlen{$readid}; $i++)
				{
					pop(@{$quals{$readid}});
					pop(@{$nqsquals{$readid}});
				}
			}
	
			### IF NEGATIVE STRAND, REVERSE ORDER OF QUAL SCORES
			if($strand{$readid} eq "-")
			{
				@{$quals{$readid}} = reverse(@{$quals{$readid}});
				@{$nqsquals{$readid}} = reverse(@{$nqsquals{$readid}});
			}
		}elsif($line =~ /[ATGC]/ && $state == 1)
		{
			$refstring .= $line;
			$state = 2;
		}elsif($line =~ /[ATGC]/ && $state == 2)
		{
			$readstring .= $line;
			$state = 1;
		}
	}
	evalRead($readstring, $refstring, $readid);
}

sub evalRead
{
	my $readstring = shift;
	my $refstring = shift;
	my $readid = shift;

	my $tempreadstring = $readstring;
	
	my $nqsqualstr = '';
	my $qualstr = '';
	
	my $lastnqs = "";
	while($tempreadstring)
	{
		my $curread = substr($tempreadstring, 0, 1);
		if($curread eq "-")
		{
			if($lastnqs == 0)
			{
				$nqsqualstr .= "9";
			}else{
				$nqsqualstr .= "8";
			}
			$qualstr .= ' ';
		}else
		{
			$qualstr .= chr(33 + shift(@{$quals{$readid}}));
			$lastnqs = shift(@{$nqsquals{$readid}});
			$nqsqualstr .= $lastnqs;
		}
		$tempreadstring = substr($tempreadstring, 1);
	}
	print QLXOUT ">Read $readid $readstart $readend ".$readlen{$readid}." ".$strand{$readid}." $contig $refstart $lastpos\n";
	print QLXOUT  "$refstring\n";
	print QLXOUT  "$readstring\n";
	print QLXOUT  "$nqsqualstr\n";
	print QLXOUT  "$qualstr\n\n";
}	
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
my $scriptpath = "/home/gbotha/galaxy-dist/tools/custom_tools/rc454";
