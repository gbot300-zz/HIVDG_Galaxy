Getting started
------------------

This package should contain the scripts runMosaik.pl, rc454.pl, qlxToSam.pl, samToQlx.pl, configPaths.pl 
and convertPos.pl.

For details on how to run any of these scripts and their options, please refer to the document 
RC454Documentation.docx in the package.

The first thing to do to get started is to run the script configPaths.pl

First, open the file configfile.txt and specify the paths where you will locate the scripts, and where
you have muscle, R (version 2.9 or higher), perl and mosaik located.

Then run the script configPaths.pl:

perl configPaths.pl configfile.txt

Testing Data
-------------------

A TestData folder is included to test the scripts and see if everything works as intended.

If you cd into this folder and run the command lines contained in VTest_Commandlines.txt, you should
in theory obtain the exact same results as the files located in the ExpectedResults folder. If you moved
the folder or the scripts in a new location make sure to modify the command lines appropriately.

Contact information
-----------------------------------------------------------------------

rc454.pl was developped from 2009 to 2011 by Patrick Charlebois for the Broad Institute of 
the MIT and Harvard

For any question or comment, contact Patrick Charlebois at patrickc@broadinstitute.org or
Matthew Henn at mhenn@broadinstitute.org
