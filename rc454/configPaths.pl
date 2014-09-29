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

my $configfile = 'configfile.txt';

open(CONFIG, $configfile);

my $scriptpath = '';
my $musclepath = '';
my $mosaikpath = '';
my $mosaik1path = '';
my $perlpath = '';
my $Rpath = '';
my $samtoolspath = '';
my $mosaiknetworkpath = '';

while(my $line = <CONFIG>)
{
	if($line =~ /scriptpath \= \'(.*)\'.*/)
	{
		$scriptpath = $1;
		if($scriptpath && substr($scriptpath, length($scriptpath) - 1) eq '/'){chop $scriptpath;}
	}elsif($line =~ /musclepath \= \'(.*)\'.*/)
	{
		$musclepath = $1;
		if($musclepath && substr($musclepath, length($musclepath) - 1) eq '/'){chop $musclepath;}
	}elsif($line =~ /mosaikpath \= \'(.*)\'.*/)
	{
		$mosaikpath = $1;
		if($mosaikpath && substr($mosaikpath && $mosaikpath, length($mosaikpath) - 1) eq '/'){chop $mosaikpath;}
	}elsif($line =~ /mosaik1path \= \'(.*)\'.*/)
	{
		$mosaik1path = $1;
		if($mosaik1path && substr($mosaik1path && $mosaik1path, length($mosaik1path) - 1) eq '/'){chop $mosaik1path;}
	}elsif($line =~ /perlpath \= \'(.*)\'.*/)
	{
		$perlpath = $1;
		if($perlpath && substr($perlpath, length($perlpath) - 1) eq '/'){chop $perlpath;}
	}elsif($line =~ /Rpath \= \'(.*)\'.*/)
	{
		$Rpath = $1;
		if($Rpath && substr($Rpath, length($Rpath) - 1) eq '/'){chop $Rpath;}
	}elsif($line =~ /samtoolspath \= \'(.*)\'.*/)
	{
		$samtoolspath = $1;
		if($samtoolspath && substr($samtoolspath, length($samtoolspath) - 1) eq '/'){chop $samtoolspath;}
	}elsif($line =~ /mosaiknetworkpath \= \'(.*)\'.*/)
	{
		$mosaiknetworkpath = $1;
		if($mosaiknetworkpath && substr($mosaiknetworkpath && $mosaiknetworkpath, length($mosaiknetworkpath) - 1) eq '/'){chop $mosaiknetworkpath;}
	}
}

open(RC454, "$scriptpath/rc454.pl");
open(RC454TMP, ">$scriptpath/rc454tmp.pl");

my $linecount = 0;
while(my $line = <RC454>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print RC454TMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print RC454TMP 'my $scriptpath = "'.$scriptpath."\";\n";
	}else{
		print RC454TMP $line;
	}
	$linecount++;
}

open(RC4541, "$scriptpath/rc454_mosaik1.pl");
open(RC454TMP1, ">$scriptpath/rc454tmp1.pl");

while(my $line = <RC4541>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print RC454TMP1 "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print RC454TMP1 'my $scriptpath = "'.$scriptpath."\";\n";
	}else{
		print RC454TMP1 $line;
	}
	$linecount++;
}


open(RUNMOSAIK, "$scriptpath/runMosaik.pl");
open(RUNMOSAIKTMP, ">$scriptpath/runMosaiktmp.pl");

$linecount = 0;
while(my $line = <RUNMOSAIK>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print RUNMOSAIKTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print RUNMOSAIKTMP 'my $scriptpath = "'.$scriptpath."\";\n";
	}elsif($line =~ /my \$mosaikpath \= .+/)
	{
		print RUNMOSAIKTMP 'my $mosaikpath = "'.$mosaik1path."\";\n";
	}else{
		print RUNMOSAIKTMP $line;
	}
	$linecount++;
}

open(RUNMOSAIK2, "$scriptpath/runMosaik2.pl");
open(RUNMOSAIKTMP2, ">$scriptpath/runMosaiktmp1.pl");

$linecount = 0;
while(my $line = <RUNMOSAIK2>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print RUNMOSAIKTMP2 "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print RUNMOSAIKTMP2 'my $scriptpath = "'.$scriptpath."\";\n";
	}elsif($line =~ /my \$mosaikpath \= .+/)
	{
		print RUNMOSAIKTMP2 'my $mosaikpath = "'.$mosaikpath."\";\n";
	}elsif($line =~ /my \$samtoolspath \= .+/)
	{
		print RUNMOSAIKTMP2 'my $samtoolspath = "'.$samtoolspath."\";\n";
	}elsif($line =~ /my \$mosaiknetworkpath \= .+/)
        {
                print RUNMOSAIKTMP2 'my $mosaiknetworkpath = "'.$mosaiknetworkpath."\";\n";
        }else{
		print RUNMOSAIKTMP2 $line;
	}
	$linecount++;
}

open(VPHASER, "$scriptpath/vphaser.pl");
open(VPHASERTMP, ">$scriptpath/vphasertmp.pl");

$linecount = 0;
while(my $line = <VPHASER>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print VPHASERTMP "#!".$perlpath." perl -w\n";
	}else{
		print VPHASERTMP $line;
	}
	$linecount++;
}

open(QLXTOSAM, "$scriptpath/qlxToSam.pl");
open(QLXTOSAMTMP, ">$scriptpath/qlxToSamtmp.pl");

$linecount = 0;
while(my $line = <QLXTOSAM>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print QLXTOSAMTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$samtoolspath \= .+/)
	{
		print QLXTOSAMTMP 'my $samtoolspath = "'.$samtoolspath."\";\n";
	}else{
		print QLXTOSAMTMP $line;
	}
	$linecount++;
}

open(SAMTOQLX, "$scriptpath/samToQlx.pl");
open(SAMTOQLXTMP, ">$scriptpath/samToQlxtmp.pl");

$linecount = 0;
while(my $line = <SAMTOQLX>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print SAMTOQLXTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$samtoolspath \= .+/)
	{
		print SAMTOQLXTMP 'my $samtoolspath = "'.$samtoolspath."\";\n";
	}else{
		print SAMTOQLXTMP $line;
	}
	$linecount++;
}


open(VPROFILER, "$scriptpath/vprofiler.pl");
open(VPROFILERTMP, ">$scriptpath/vprofilertmp.pl");

$linecount = 0;
while(my $line = <VPROFILER>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print VPROFILERTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$scriptpath \= .+/)
	{
		print VPROFILERTMP 'my $scriptpath = "'.$scriptpath."\";\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print VPROFILERTMP 'my $musclepath = "'.$musclepath."\";\n";
	}elsif($line =~ /my \$Rpath \= .+/)
	{
		print VPROFILERTMP 'my $Rpath = "'.$Rpath."\";\n";
	}else{
		print VPROFILERTMP $line;
	}
	$linecount++;
}

open(CONVERTPOS, "$scriptpath/convertPos.pl");
open(CONVERTPOSTMP, ">$scriptpath/convertPostmp.pl");

$linecount = 0;
while(my $line = <CONVERTPOS>)
{
	if($linecount == 0 && $line =~ /^\#\!.+/)
	{
		print CONVERTPOSTMP "#!".$perlpath." perl -w\n";
	}elsif($line =~ /my \$musclepath \= .+/)
	{
		print CONVERTPOSTMP 'my $musclepath = "'.$musclepath."\";\n";
	}else{
		print CONVERTPOSTMP $line;
	}
	$linecount++;
}

system("mv $scriptpath/rc454tmp.pl $scriptpath/rc454.pl");
system("mv $scriptpath/rc454tmp1.pl $scriptpath/rc454_mosaik1.pl");
system("mv $scriptpath/convertPostmp.pl $scriptpath/convertPos.pl");
system("mv $scriptpath/runMosaiktmp.pl $scriptpath/runMosaik.pl");
system("mv $scriptpath/runMosaiktmp1.pl $scriptpath/runMosaik2.pl");
system("mv $scriptpath/vphasertmp.pl $scriptpath/vphaser.pl");
system("mv $scriptpath/vprofilertmp.pl $scriptpath/vprofiler.pl");
system("mv $scriptpath/qlxToSamtmp.pl $scriptpath/qlxToSam.pl");
system("mv $scriptpath/samToQlxtmp.pl $scriptpath/samToQlx.pl");

print "Config done\n";
