#!/usr/bin/bash

PAD=25
BASE=`echo $1 | cut -f 1 -d '.'` #assumes a file *.bed
WNAME=`echo $BASE | perl -ne "chomp; print; print \"_withSeqs.txt\"; "`

echo "Processing file:$1 through base:$BASE to create countWithSeqFile:$WNAME"

~/prog/faUtil/faGetBedSequences /Users/jhgraber/data/yeast/sacCer3_genome_plus_2mic_ERCC92_one-line-seqs.fa $1 | perl -ne "s/_/\t/g; print" | perl -ane "print \"@F[0]\t@F[1]\t@F[2]\t@F[3]\t@F[4]\t@F[5]\t\"; print substr(@F[6],0,$PAD); print \"\t\"; print substr(@F[6],$PAD,2*$PAD); print \"\t@F[7]_\"; print ((@F[2]+@F[3])/2); print \"_@F[4]_@F[5]\t@F[12]\t@F[13]\n\"; " > $WNAME

