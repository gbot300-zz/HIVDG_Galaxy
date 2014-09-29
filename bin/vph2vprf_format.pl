#========================================================================
# Project     : ExomeHLATyper
# Name        : vph2vprf_format.pl
# Author      : Xiao Yang
# Created on  : March, 2013
# Version     : 1.0
# Copyright   : The Broad Institute
#  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
# 				 This software and its documentation are copyright (2012)
#	  			 by the Broad Institute. All rights are reserved.
#
# 				 This software is supplied without any warranty or 
#		  		 guaranteed support whatsoever. The Broad Institute cannot 
#				   be responsible for its use,	misuse, or functionality.
# Description :
#========================================================================


use strict;
use Getopt::Long;
use Data::Dumper;
 
my %option = (
	h           => '',
	ref					=> '',
	vph2				=> '',
	o						=> '',
);

GetOptions(
	"h"						=> \$option{h},
	"ref=s"				=> \$option{ref},
	"vph2=s"			=> \$option{vph2},
	"o=s"					=> \$option{o},	
) || printHelp (); 

if ($option{h}) { printHelp();}


unless ($option{ref} && $option{vph2} && $option{o}) {	printHelp(); }

sub printHelp {
		print "\n--------------------------------------------------------------------\n";
		print "usage: perl vph2vprf_format.pl -ref [ref.fasta] -vph2 [vphaser2.output.txt] -o [vprofiler_format.txt]\n\n";
		
		print "\t-ref: reference genome\n";
		print "\t-vph2: vphaser 2.0 output variant file\n";				
		print "\t-o: vprofiler input format\n";				
		print "--------------------------------------------------------------------\n\n";
		exit;
}
 
# read in reference genome
open (REF, "<$option{ref}") or die "unable to open file $option{ref} to read\n";
my $seq;
while (<REF>){
	if (/>/) {
	} else {
		s/\s+//g; # remove whitespace
    $seq .= $_; # add sequence
	}
}
close (REF);
uc($seq); # convert to upper case

# read in vphaser 2.0 output variant file
my %variants;

open (VPH, "<$option{vph2}") or die "unable to open file $option{vph2} to read\n";

while (<VPH>) {
	my $line = $_;
	if ($line !~ /#/) {
		my @entries = split ("\t",$line);
		my $num = $#entries + 1;
		if ($num < 3) { die "line - $line - in $option{vph2} has less than 3 entries" };
		my $pos = $entries[0];
		my $var = $entries[1];
		
		if ($var eq "i") { # change var to cons, this call is wrt to the ref genome
			$var = $entries[2];
		}
		
#		print $var."\n";
#		if ($var eq "C") { print "ok";}
#		exit;
		my $cons = $entries[2];

		# check if pos is already registered
		if (!$variants{$pos}) {
			my @tmp = (0, 0, 0, 0, 0, 0);
			$variants{$pos} = \@tmp ;
		} 
#		print "variants = $variants{$pos}->[1]\n";
#		exit;

		if ($var eq "A" || $cons eq "A") { $variants{$pos}->[0] = 1;	} 
		if ($var eq "T" || $cons eq "T") { $variants{$pos}->[1] = 1;	} 
		if ($var eq "G" || $cons eq "G") { $variants{$pos}->[2] = 1;	} 
		if ($var eq 'C' || $cons eq "C") { $variants{$pos}->[3] = 1;	}
		
		my $num_del = 0;
		if ($var =~ /D(\d+)/) {	
			$num_del = $1;
			if ($variants{$pos}->[4] > 0) { # more types of deletions found, push to the end of array
				push (@{$variants{$pos}}, $num_del);  
			} else {
				$variants{$pos}->[4] = $num_del;	
			}
#			-- $num_del;
		}	
		
		if ($var =~ /I/) {	$variants{$pos}->[5] = 1;	}
		if ($var =~ /I/ || $var =~ /D/) { # identify the corresponding ref base
			my $base = substr($seq, $pos - 1, 1);
			my $index = 0; # default A for consensus 
			if ($base eq "T") { $index = 1; }
			elsif ($base eq "G") { $index =	2; }
			elsif ($base eq "C") { $index = 3; }
			$variants{$pos}->[$index] = 1;			
			
#			if ($num_del > 0) {	
#				for (my $i = 1; $i <= $num_del; ++ $i) { # identify how long this deletion is 
#					my @tmp = (0, 0, 0, 0, 0, 0);
#					$variants{$pos + $i} = \@tmp ;			
#					$variants{$pos + $i}->[4] = 1;
#					$variants{$pos + $i}->[$index] = 1;
#				}
#			}
			
		} # if ($var =~ /I/ || $var =~ /D/)
	} # ($line !~ /#/)
} # while (<VPH>)

close (VPH);

# write output
open (OUTPUT, ">$option{o}") or die "unable to open file $option{o} to write\n";
print OUTPUT "Position\tA\tT\tG\tC\tDeletion\tInsertion\n";

my $ref_len = length ($seq);

for (my $index = 0; $index < $ref_len; ++ $index) {
	my $i = $index + 1;
	if (!$variants{$i}) { # variant doesn't exists
		my $base = substr($seq, $index, 1);
		if ($base eq "A") {
				print OUTPUT "$i\t1\t0\t0\t0\t0\t0\n";
		} elsif ($base eq "T") {
				print OUTPUT "$i\t0\t1\t0\t0\t0\t0\n";		
		} elsif ($base eq "G") {
				print OUTPUT "$i\t0\t0\t1\t0\t0\t0\n";		
		} elsif ($base eq "C") {
				print OUTPUT "$i\t0\t0\t0\t1\t0\t0\n";
		}
	} else { # variant
		my @deletions;
		my $variant_num =  scalar (@{$variants{$i}}) ;
#		print $variant_num."\n";
#		exit;
		if ($variant_num > 6) {
			for (my $j = $variant_num - 1; $j > 5; ++ $j) {
				push (@deletions, $variants{$i}->[$j]);
			}	
		} 
		
#		print OUTPUT "$i\t$variants{$i}->[0]\t$variants{$i}->[1]\t$variants{$i}->[2]\t$variants{$i}->[3]\t$variants{$i}->[4]\t$variants{$i}->[5]\t\n";
		
		print OUTPUT "$i\t$variants{$i}->[0]\t$variants{$i}->[1]\t$variants{$i}->[2]\t$variants{$i}->[3]\t$variants{$i}->[4]\t";
		foreach (@deletions) {
 				print OUTPUT $_."\t";
 		} 
 		print OUTPUT "$variants{$i}->[5]\n";
	}
}
close (OUTPUT);
