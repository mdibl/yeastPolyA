#!/usr/bin/bash

PAD=25
BASE=`echo $1 | cut -f 1 -d '.'`
BNAME=`echo $BASE | perl -ne "chomp; print; print \".bed\"; "`
SNAME=`echo $BASE | cut -f 1 -d '_'`
echo "Processing file:$1 through base:$BASE to create bedfile:$BNAME, with sample:$SNAME"

perl -ne "s/_\+/_#/; s/_\-/_+/; s/_#/_-/; s/_/\t/g; print; " $1 | perl -ane "print \"@F[0]\t\"; print (@F[1]-$PAD); print \"\t\"; print (@F[1]+$PAD); print \"\t$SNAME\t@F[3]\t@F[2]\t@F[4]\t@F[5]\n\"; " > $BNAME
