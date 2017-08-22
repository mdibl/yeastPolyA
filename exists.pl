#!/usr/bin/perl

# Script Name: exists.pl
# Author: Kristoph S. N. Naggert, Wednesday August 9th 2017
# Type: Perl Script

# Arguments: Script passes in a design file as the only argument.

# Execution: perl exists.pl (designfile) 

# Overview: This script is used to read in a specific design file containing the information for the project in the form: prefix (prefix of .psl file name), label, mut (mutation i.e. wild type, ipa1-1), batch, rep; it then checks for the existence of the .psl files with the prefix read in from the design file. If those files do not exists the program quits and prints and error message informing the user that those files do not exist in the specified directory, if the files do exists the script creates a copy of the base makefile and appends several rules to the end of it needed to create the datafile used in the final analysis. 

# Variables: There are a number of variables the will need to be replaced by the user, ideally this script will be updated at some point to move those variables into a configuration file that can be manipulated by user of the program.

# $ARGV[0] = Design file for the project (Note: This will not be defined in the script but rather as the only argument you pass in at the command line)

# my $psllocation = User directory where .psl files are located

# my $makefilelocation = User directory where the base makefile is located

# my $makefileBase = Name of base makefile (Note: All makefiles must be named ‘makefile’ in order to be executed, by naming the base set of rules something other than that will have the cpu recognize it as a simple text file, either type of files should work for this, however this was tested in a setting where the base rules were recognized as a text file.)

# my $assignedMakefile = Name of the makefile you are creating (Note: This needs to be ‘makefile’ otherwise the cpu will not recognize the file as a makefile, and it will not be able to be run as such.)

# my $existingdir = User directory where the new makefile will be created (Note: Its easiest if this is the same directory where the files created by the makefile are being saved to.)


use strict;
use warnings;
use FileHandle;
use File::Copy;


my $psllocation = "/Users/knaggert/Desktop/pcf11_summer/abc/";
my $makefilelocation = "/Users/knaggert/Desktop/pcf11_summer/";
my $makefileBase = "makefile_base";
my $assignedMakefile = "makefile";
my $existingdir = "/Users/knaggert/Desktop/pcf11_summer/outFiles";

my $fh = FileHandle->new;

print $ARGV[0], "\n"; #This should print the name of the design file
open($fh, $ARGV[0]) || die "Cannot open $ARGV[0]: $!"; #Attempts to Open File, if it Can’t Kills Program and Returns Error Message

while(my $line = <$fh>) { #While There are Still Lines in the File Loop Through Them One at a Time
	my @prefix = split /\s+/,$line; #Split Up Each Line into Individual Words [Looks for One or More White Spaces Just in Case the Design File is not Tab Delimitated]
	if ($prefix[0] eq "prefix") { #Check if you are in the Header Row
		print "This is the header row, no need to check for a .psl file here.\n\n";
		}
	else {
		my $pslFile = $prefix[0].'.psl'; #Takes the Prefix[First Entry in the @Prefix Variable] from the Design File and Adds the .psl Tag
		if (-e "$psllocation/$pslFile") {
			print "File $psllocation/$pslFile exists.\n\n"; #Prints out the location of the file if it exists
			}
		else {
			die "$pslFile does not exist in $psllocation"; #If any of the .psl files do not exist the the program is killed and an error message is returned
			}
		
	}
}
close($fh);

#This builds the $assignedFileRule which will be appended to the makefile

my $assignedFileRule = "";

my $fileH = FileHandle->new;

open($fileH, $ARGV[0]) || die "Cannot open $ARGV[0]: $!"; #Attempts to open the design file, if it can't it kills the program and returns an error message

while (my $line = <$fileH>) { #While There are Still Lines in the File Loop Through Them One at a Time
	my @prefix = split /\s+/, $line; #Split Up Each Line into Individual Words [Looks for One or More White Spaces Just in Case the Design File is not Tab Delimitated]
	my $assignedFile = $prefix[0].'_assigned.txt'; #sets the first element in the row equal to the assignedFile variable
	if ($prefix[0] eq "prefix") { #Check if you are in the Header Row
		}
	else {
		if (grep{$_ eq $assignedFile} $assignedFileRule) { #If the file has already been added to the rule it returns a message saying that
			print "This file is already in the list";
		}
		else { #If the file is not already in the rule this adds it to that list
			$assignedFileRule = $assignedFileRule." ".$assignedFile;
		}
	}
}

print "$assignedFileRule";


my $fileHandle = FileHandle->new;

if (-e "$makefilelocation/$makefileBase") { #Checks to see if base makefile exists in the specified directory
	print "\n\nFile $makefilelocation/$makefileBase exists.\n"; #Prints out the location of the base makefile if it exists
	copy("$makefileBase", "$existingdir/$assignedMakefile") or die "Copy failed: $!"; #Creates a copy of the existing makefile, if it cannot the program is killed and an error message is returned
	open (my $fileHandle, '>>', "$existingdir/$assignedMakefile") or die "Can't open '$existingdir/assignedMakefile'\n"; #Opens the new makefile, if it can't it kills the program and returns an error message
	#Appends the new set of rules to the makefile
	print $fileHandle "\n\nassignedCat.txt: $assignedFileRule\n";
	print $fileHandle "\tcat \$^ > \$(expNm)_\$@\n\n";
	print $fileHandle "CollectSiteCountsOutputFiles: \$(expNm)_assignedCat.txt \$(CPATH)CollectSiteCounts\n";
	print $fileHandle "\t\$(CPATH)CollectSiteCounts \$< \$(expNm)\n\n";
	print $fileHandle "processAllPolyAFiles: \$(expNm)_site.txt \$(CPATH)processAllPolyASites\n";
	print $fileHandle "\t\$(CPATH)processAllPolyASites \$< \$(expNm)\n\n";
	close $fileHandle;
			}
else { #If the File Can’t be Found it Kills the Program and Returns an Error Message
	die "$makefileBase does not exist in $makefilelocation";
		}
