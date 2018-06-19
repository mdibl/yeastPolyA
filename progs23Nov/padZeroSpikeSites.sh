#!/usr/bin/bash

#!/usr/bin/bash

BASE=`echo $1 | cut -f 1 -d '.'`
echo Processing $BASE
ZNAME=`echo $BASE | perl -ne "chomp; print; print \"_wZeros.txt\"; "`

join -v 1 ERCC_allSites.txt $1 | perl -ne "chomp; print; print \"\t0\t0\t0\n\"; " | cat - $1 | sort > $ZNAME
 