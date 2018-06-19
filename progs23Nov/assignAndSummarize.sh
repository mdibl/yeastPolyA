#!/usr/bin/bash

BASE=`echo $1 | perl -ne "chomp; s/_siteCount.*//; print;"` 
ANAME=`echo $BASE | perl -ne "chomp; print; print \"_assigned.txt\"; "`

echo "Processing file:$1 through base:$BASE to create assignedFile:$ANAME"
/Users/jhgraber/data/yeast/progs23Nov/AssignSitesToGenes $1 > $ANAME
echo "Summarizing assginedFile:$ANAME to base:$BASE"
/Users/jhgraber/data/yeast/progs23Nov/SummarizeGeneCounts $ANAME $BASE

