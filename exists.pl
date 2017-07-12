#Perl Script designed to read in a specific design file containing samples and checks for the existence of psl files for all samples.
#!/usr/bin/perl

#Kristoph Naggert, Thursday July 6th, 2017

#Loop through each line of the design.txt file [While Loop]
#In each line split the line and take the first word in that line(Prefix of the .psl file) [Split]
#Take the first word(Prefix of the .psl file) and add the .psl tag [Join]
#Conditional statement, if the .psl file exists, continue else exit the program and return error [If/Else]



use strict;
use warnings;

$filename = $ARGV[0]; #Tab Delimited Design File
$psllocation = /Users/knaggert/Desktop/pcf11_summer/yeastData/; #User Directory Where .psl Files are Kept
$makefilelocation = /Users/knaggert/Desktop/pcf11_summer/;
$makefileBase = $ARGV[1]; #Base Makefile
$existingdir = /Users/knaggert/Desktop/pcf11_summer/outFiles’; #Directory where the makefile will be created

open($fh => $filename) || die "Cannot open $filename: $!"; #Attempts to Open File, if it Can’t Kills Program and Returns Error Message

while($line = <$fh>) { #While There are Still Lines in the File Loop Through Them One at a Time
	if (@prefix[0] == “prefix”) { #Check if you are in the Header Row
		print “This is the header row, no need to check for a .psl file here”;
	}
	else {
    	@prefix = split(“/\s+/”,$line); #Split Up Each Line into Individual Words [Looks for One or More White Spaces Just in Case the Design File is not Tab Delimitated]
		$pslFile = join(“.psl”, @prefix[0]); #Takes the Prefix[First Entry in the @Prefix Variable] from the Design File and Adds the .psl Tag
		if (-e "$location/$pslFile”) {
    		print "File $location/$file exists.\n";
			}
		else { 
		     die “$pslFile does not exist in $location”;
		}
	}
}
}
close($fh);

if (-e "$location/$pslFile”) { #If the File Exists Prints That the File Exists
    		print "File $makefilelocation/$makefileBase exists.\n";
		open my $fileHandle, ">>", "$existingdir/assignedMakefile” or die "Can't open '$existingdir/assignedMakefile'\n"; 
		system(cat makefileBase >> assignedMakefile);
		print “assignedCat.txt : $(assignedFiles)\n”;
		print ““\t”cat $(assignedFiles) > $@\n\n”;
		
		close $fileHandle;
			}
else { #If the File Can’t be Found it Kills the Program and Returns an Error Message
			die “$makefileBase does not exist in $makefilelocation”;
		}


