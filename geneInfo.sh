#!/usr/bin/bash

# Script Name: geneInfo.sh
# Author: Kristoph S. N. Naggert, Wednesday August 9th 2017
# Type: BASH Script

# Arguments: Script passes in a gene name, a cumPa.txt file, a feature file, a Chip-seq file, and a summarized feature file

# Execution: bash geneInfo.sh (gene name) (cumPa.txt file) (feature file) (Chip-seq file) (summarized feature file)

# Overview: This script gathers all of the pieces of information needed to create plots of cumulative poly(A) probability, Pol II enrichment, and difference in Pol II enrichment for a specific gene and puts them into specific files. It then builds the gnuplot script to build the graphics and then runs those scripts.

# Variables: There are a number of variables the will need to be replaced by the user.

# $1 = The name of the gene you are looking for (passed in at the command line)

# $2 = The name of the cumPa.txt file being searched (passed in at the command line)

# $3 = The name of the feature file being searched (passed in at the command line)

# $4 = The name of the chip-Seq file being searched (passed in at the command line)

# $5 = The summarized feature file being searched (passed in at the command line)

# upstreamBuffer = The number of base pairs upstream of the gene.

# downstreamBuffer = The number of base pairs upstream of the gene; Ideally both downstreamBuffer and upstreamBuffer will be moved into a configuration file that can be manipulated from outside the script.

# start = The start position of the gene passed in as the first argument (Found in the 3rd column of the feature file)

# stop = The stop position of the gene passed in as the first arguement (Found in the 4th column of the feature file)

# SEN = The direction of the gene passed in as the first argument (Found in the 5th column of the feature file)

# right_chromosome = The chromosome of the gene passed in as the first argument (Found in the 1st column of the feature file)

# offset = The difference in Pol II enrichment between Wild Type and Mutant (Found in the 11th column of the summarized feature file)

printf "The gene being searched for is $1\nThe *_cumPa.txt file being searched is $2\nThe feature file being searched is $3\nThe Chip-seq file being searched is $4\nThe featureSummarized file being searched is $5\n\n";

upstreamBuffer=100;
downstreamBuffer=500;

grep $1 $2 > $1_cumPa.txt;

start=`grep $1 $3 | cut -f 3`
echo $start
stop=`grep $1 $3 | cut -f 4`
echo $stop
SEN=`grep $1 $3 | cut -f 5`
right_chromosome=`grep $1 $3 | cut -f 1`
offset=`grep $1 $5 | cut -f 11`

#Checks the orientation of the gene and then applies the up and downstream buffers appropriatly
if [ $SEN == "+" ];  then
	start=`echo $start | perl -ane "print(@F[0]-$upstreamBuffer);"`
	stop=`echo $stop | perl -ane "print(@F[0]+$downstreamBuffer);"`
else
	start=`echo $start | perl -ane "print(@F[0]-$downstreamBuffer);"`
	stop=`echo $stop | perl -ane "print(@F[0]+$upstreamBuffer);"`
fi

#Appends the Pol II enrichment 'offset' to the gene specific chipSeq file.
perl -ane "if(/start/) {chomp; print; print \"\t\"; print 'offset'; print \"\n\"; } elsif(@F[0] eq $right_chromosome && @F[1] >= $start && @F[1] <= $stop) {chomp; print; print \"\t\"; print (@F[9]- $offset); print \"\n\"; }" $4 > $1_chipSeq.txt;

#Checks the orientation of the gene, if the gene is antisense then if sets the REV variable to "reverse" so in the plots it will be plotted in reverse.
if [ $SEN == "-" ]; then
	echo plotting antisense
	REV="reverse"
fi
#This block builds the GNUplot script that will creates the plots we've specified. It starts by setting the .png size, and then it sets the name of the output png, genename_plots.png. Next it sets up multiplot which allows us to combine several plots into one. Next we set up the key for the first plot. Then the title is set for the first plot, the y axis is labeled, the y and x tics fonts are set to Arial narrow, and xrange is reversed because REV = reverse. And finally we plot, using the genename_cumPa.txt file. Then we unset the title, the yaxis label and the key. For the next plot we create another key, the title is set, the y axis label is set and then we plot using the genename_chipSeq.txt file. Then we unset the title, ylabel and key. For the final plot we once again set the key, then we set the title for the third plot, the we set the y axis label and plot using the genename_chipSeq.txt file. We unset multiplot and save the gnuplot script to a file called genename_script.txt, then we pipe that script into gnuplot and open genename_plot.png.
echo "set term png size 1200,800 font \"/Library/Fonts/Arial Narrow.tff,10\"
set output '$1_plots.png'
set multiplot layout 3, 1 title '$1' font '/Library/Fonts/Arial Narrows.ttf,16'
set key outside top right vertical Right width 0 spacing 2 box font '/Library/Fonts/Arial Narrows.ttf,9'
set title '$1 Cumulative Poly(A)' offset 0,-1 font '/Library/Fonts/Arial Narrow.ttf,14'
set ylabel 'cumulative poly(A) probability' font '/Library/Fonts/Arial Narrow.ttf,14'
set xtics font '/Library/Fonts/Arial Narrow.ttf,11' 
set ytics font '/Library/Fonts/Arial Narrow.ttf,11' 
set xrange [] $REV writeback
set rmargin 30
plot '$1_cumPa.txt' u 3:11 t 'Lane8ByIB' w lines lc rgb '#dda0dd' lw 1, '$1_cumPa.txt' u 3:12 t 'Lane8ByIIB' w lines lc rgb '#dda0dd' lw 1, '$1_cumPa.txt' u 3:13 t 'Lane8ByIIIB' w lines lc rgb '#dda0dd' lw 1, '$1_cumPa.txt' u 3:14 t 'Lane8TS1248IB' w lines lc rgb '#ffa500' lw 1, '$1_cumPa.txt' u 3:15 t 'Lane8TS1248IIIB' w lines lc rgb '#ffa500' lw 1, '$1_cumPa.txt' u 3:16 t 'lane3BY4741III' w lines lc rgb '#dda0dd' lw 1, '$1_cumPa.txt' u 3:17 t 'lane3TS1248II' w lines lc rgb '#ffa500' lw 1, '$1_cumPa.txt' u 3:(0.25*(\$11+\$12+\$13+\$16)) t 'WildType Average' w lines lc rgb '#9400d3' lw 2, '$1_cumPa.txt' u 3:(0.33*(\$14+\$15+\$17)) t 'ipa1-1 Average' w lines lc rgb '#ff4500' lw 2
unset title
unset ylabel
unset key
set xrange restore
set key outside top right vertical Right width 0 spacing 2 box font '/Library/Fonts/Arial Narrows.ttf,9'
set title '$1 Pol-II Enrichment' offset 0,-1 font '/Library/Fonts/Arial Narrow.ttf,14'
set ylabel 'Pol-II Enrichment' font '/Library/Fonts/Arial Narrow.ttf,14'
plot '$1_chipSeq.txt' u 2:3 t 'WT_log2_enrichment_rep2' w lines lc rgb '#bebebe' lw 2, '$1_chipSeq.txt' u 2:4 t 'WT_log2_enrichment_rep1' w lines lc rgb '#bebebe' lw 2, '$1_chipSeq.txt' u 2:5 t 'Ipa1-1_log2_enrichment_rep2' w lines lc rgb '#ffc0c0' lw 2, '$1_chipSeq.txt' u 2:6 t 'Ipa1-1_log2_enrichment_rep1' w lines lc rgb '#ffc0c0' lw 2, '$1_chipSeq.txt' u 2:7 t 'WT' w lines lc rgb '#000000' lw 2, '$1_chipSeq.txt' u 2:8 t 'ipa1-1' w lines lc rgb '#ff0000' lw 2
unset title
unset ylabel
unset key
set key outside top right vertical Right width 0 spacing 2 box font '/Library/Fonts/Arial Narrows.ttf,9'
set title '$1 Difference in Pol-II Enrichment' offset 0,-1 font '/Library/Fonts/Arial Narrow.ttf,14'
set ylabel 'Difference in Pol-II Enrichment' font '/Library/Fonts/Arial Narrow.ttf,14'
plot '$1_chipSeq.txt' u 2:9 t 'Difference' w line lc rgb '#ee82ee' lw 2, '$1_chipSeq.txt' u 2:(\$9+\$10) t 'Difference with Offset' w line lc rgb '#c080ff' lw 2
unset multiplot" > $1_script.txt
cat $1_script.txt | gnuplot
open $1_plots.png

