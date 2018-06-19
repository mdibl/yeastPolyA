#!/usr/bin/bash

BASE=`echo $1 | cut -f 1 -d '.'`
echo Processing $BASE
ANAME=`echo $BASE | perl -ne "chomp; print; print \"_alnCount.txt\"; "`
TNAME=`echo $BASE | perl -ne "chomp; print; print \"_tgtCount.txt\"; "`
JNAME=`echo $BASE | perl -ne "chomp; print; print \"_tgt_aln.txt\"; "`
CNAME=`echo $BASE | perl -ne "chomp; print; print \"_siteCount.txt\"; "`

grep seq $1 | perl -ane "if(@F[0]>=16 && @F[17]<2) { print \"@F[9]\n\";} " | uniq -c | perl -ane "print \"@F[1]\t@F[0]\n\"; " | sort > $ANAME
grep seq $1 | perl -ane "if(@F[0]>=16 && @F[17]<2) { if(@F[8] eq \"-\") { print \"@F[13]_@F[16]_-\n\";} else { print \"@F[13]_@F[15]_+\n\";} } " | sort | uniq -c | perl -ane "print \"@F[1]\t@F[0]\n\"; " > $TNAME
grep seq $1 | perl -ane "if(@F[0]>=16 && @F[17]<2) { if(@F[8] eq \"-\") { print \"@F[9]\t@F[13]_@F[16]_-\n\";} else { print \"@F[9]\t@F[13]_@F[15]_+\n\";} } " | sort -k 2 | join -1 2 - $TNAME | perl -ane "print \"@F[1]\t@F[0]\t@F[2]\n\"; " | sort | join - $ANAME | perl -ane "print \"@F[1]\t@F[0]\t@F[2]\t@F[3]\n\"" > $JNAME
~/prog/yeast-pA-flow/CondenseTagsToSiteCount $JNAME > $CNAME
