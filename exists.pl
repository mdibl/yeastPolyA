#Perl Script designed to read in a specific design file containing samples and checks for the existence of psl files for all samples.
#!/usr/bin/perl

#Kristoph Naggert, Thursday July 6th, 2017

#Loop through each line of the design.txt file [While Loop]
#In each line split the line and take the first word in that line(Prefix of the .psl file) [Split]
#Take the first word(Prefix of the .psl file) and add the .psl tag [Join]
#Conditional statement, if the .psl file exists, continue else exit the program and return error [If/Else]



use strict;
use warnings;
use FileHandle;
use File::Copy;

#my $filename = pcf11_caff_design.txt; #Tab Delimited Design File
my $psllocation = "/Users/knaggert/Desktop/pcf11_summer/yeastData/"; #User Directory Where .psl Files are Kept
my $makefilelocation = "/Users/knaggert/Desktop/pcf11_summer/";
my $makefileBase = "makefile_base"; #Base Makefile
my $assignedMakefile = "makefile";
my $existingdir = "/Users/knaggert/Desktop/pcf11_summer/outFiles"; #Directory where the makefile will be created

my $fh = FileHandle->new;

print $ARGV[0], "\n"; 
open($fh, $ARGV[0]) || die "Cannot open $ARGV[0]: $!"; #Attempts to Open File, if it Can’t Kills Program and Returns Error Message

while(my $line = <$fh>) { #While There are Still Lines in the File Loop Through Them One at a Time
	my @prefix = split /\s+/,$line; #Split Up Each Line into Individual Words [Looks for One or More White Spaces Just in Case the Design File is not Tab Delimitated]
	if ($prefix[0] eq "prefix") { #Check if you are in the Header Row
		print "This is the header row, no need to check for a .psl file here";
		}
	else {
		my $pslFile = $prefix[0].'.psl'; #Takes the Prefix[First Entry in the @Prefix Variable] from the Design File and Adds the .psl Tag
		if (-e "$psllocation/$pslFile") {
			print "File $psllocation/$pslFile exists.\n\n";
			}
		else {
			die "$pslFile does not exist in $psllocation";
			}
		
	}
}
close($fh);

my $assignedFileRule = "";

my $fileH = FileHandle->new;

open($fileH, $ARGV[0]) || die "Cannot open $ARGV[0]: $!";

while (my $line = <$fileH>) {
	my @prefix = split /\s+/, $line;
	my $assignedFile = $prefix[0].'_assigned.txt';
	if ($prefix[0] eq "prefix") { #Check if you are in the Header Row
		}
	else {
		if (grep{$_ eq $assignedFile} $assignedFileRule) {
			print "This file is already in the list";
		}
		else {
			$assignedFileRule = $assignedFileRule." ".$assignedFile;
		}
	}
}

print "$assignedFileRule";


my $fileHandle = FileHandle->new;

if (-e "$makefilelocation/$makefileBase") {
	print "\n\nFile $makefilelocation/$makefileBase exists.\n";
	copy("$makefileBase", "$existingdir/$assignedMakefile") or die "Copy failed: $!";
	open (my $fileHandle, '>>', "$existingdir/$assignedMakefile") or die "Can't open '$existingdir/assignedMakefile'\n";
	print $fileHandle "assignedCat.txt: $assignedFileRule\n"; #Just prints to console, incorrect syntax
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
