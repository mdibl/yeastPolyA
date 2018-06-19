#!/usr/bin/bash

FILES=*siteCount.txt
for f in $FILES 
do
    BASE=`echo $f | cut -f 1 -d '.'`
    ZNAME=`echo $BASE | perl -ne "chomp; print; print \"_wZeros.txt\";"`
    #ONAME=`echo $BASE | perl -ne "chomp; print; print \"_oneCount.txt\";"`
    echo Processing $f

    join -v 1 genome_allSites_atLeastOneSampGT0.txt $f | perl -ne "chomp; print; print \"\t0\t0\t0\n\"; " | cat - $f | sort > $ZNAME
done