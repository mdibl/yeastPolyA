grep -E "$1$"  ../saccer3_simple_features.txt | cut -f 3,4 | perl -ane "print \"@F[0]\t-0.03\n@F[1]\t-0.03\n\";" > plots/$1_coord.txt 
grep -E "$1[[:space:]]" 4nqo_std_Arem_unique_1000_05.edges.txt | grep -v totals | grep -v bounds | grep -v STD > plots/$1.txt
#grep tota plots/$1.txt | cut -f 3-4 | perl -ane "print \"STD\t\"; print (1*@F[1]/476); print \"\n4NQO\t\"; print (1*@F[0]/1363); print \"\n\"; " > plots/$1_exp.txt
grep -E "$1[[:space:]]" 4nqoExpComp.txt | cut -f 2,5 | perl -ane "print \"STD\t@F[0]\n4NQO\t@F[1]\n\"; " > plots/$1_exp_tot.txt
grep -E "$1[[:space:]]" 4nqoExpComp.txt | cut -f 3,6 | perl -ane "print \"STD\t@F[0]\n4NQO\t@F[1]\n\"; " > plots/$1_exp_fl.txt
grep -E "$1[[:space:]]" 4nqoExpComp.txt | cut -f 4,7 | perl -ane "print \"STD\t@F[0]\n4NQO\t@F[1]\n\"; " > plots/$1_exp_use.txt
grep -E "$1[[:space:]]" 4nqo_std_Arem_unique_30_10.out.txt | cut -f 7,8 | perl -ane "print \"@F[0]\t-0.06\n@F[1]\t-0.06\n\";" > plots/$1_block.txt
#./efTest plots/$1.txt 1 3 12 0.05 100 > plots/$1_edge.txt
#echo "set term aqua font 'Arial,16' size 900 600 title '$1' 
PVAL=`grep -E "$1[[:space:]]" 4nqo_std_Arem_unique_30_10.out.txt | cut -f 5` 
PRT=`grep -E "$1[[:space:]]" 4nqo_std_Arem_unique_30_10.out.txt | cut -f 5 | perl -ane "if(@F[0]<0.05) { print \"yes\"; } else { print \"no\"; }"`
SEN=`grep -E "$1$"  ../saccer3_simple_features.txt | cut -f 5`
SIG=`cut -f 8 plots/$1.txt | uniq | wc -l`
REV="noreverse"
BLK=""
SGP=""
if [ $SEN == "-" ]; then
   echo plotting antisense
   REV="reverse"   
fi
if [ $PRT == "yes" ]; then
   echo block found
   BLK=", 'plots/$1_block.txt' using 1:2 with l lc rgb '#9999FF' lw 6 t 'MaxBlock'"
fi
if [ $SIG != "1" ]; then
   echo sig sites found
   SGP=", 'plots/$1.txt' using 2:(\$8-1.08) w points pt 1 ps 1 lc rgb 'grey' t 'SigBlocks'"
fi
echo "set term png size 1000 600 font \"/Library/Fonts/Arial Rounded Bold.ttf,10\"  
set output 'plots/$1_4nqo.png'
#set key outside bottom center horizontal box
set key inside top center horizontal box  font \"/Library/Fonts/Arial Rounded Bold.ttf,9\" 
set multiplot
set size 0.9, 0.75
set origin 0.0, 0.25
set grid lt 0
set title 'Raw Sites $1: `grep -E $1[[:space:]] alias_commas.txt | cut -f 2`' offset 0,-0.5 font '/Library/Fonts/Arial Rounded Bold.ttf,16'
set ylabel 'cumulative poly(A) probability' offset 3 font '/Library/Fonts/Arial Rounded Bold.ttf,11' 
set ytics font '/Library/Fonts/Arial Rounded Bold.ttf,10' offset 1
set bmargin 0
set xtics format ''
set yrange [-0.1:1.1]
set xrange [] $REV writeback
set rmargin 15
  plot \
     \"plots/$1.txt\" using 2:4 with l lc rgb 'black' lw 2 t 'STD',\
     \"plots/$1.txt\" using 2:3 with l lc rgb 'red' lw 2 t '4NQO',\
     \"plots/$1_coord.txt\" using 1:2 with l lt 9 lw 8 t 'CDS' $BLK $SGP
set xtics format \"%.0f\" font '/Library/Fonts/Arial Rounded Bold.ttf,10' offset 0,0.5
set xlabel 'position (nt)' font '/Library/Fonts/Arial Rounded Bold.ttf,11' offset 0,1
set size 0.9,0.245
set origin 0.0,0.0
unset bmargin
unset title
set ylabel 'difference' font '/Library/Fonts/Arial Rounded Bold.ttf,11'
set yrange [-0.6:0.6]
set y2tics format ''
set y2range [-0.1:0.1]
set ytics -0.5,0.5
set nokey
set xrange restore
set tmargin 0.2
plot \"plots/$1.txt\" using 2:6 with filledc above y1=0 lc rgb '#CCAAAA', '' using 2:6 with filledc below y1=0 lc rgb '#AAAACC', '' using 2:5 with l lw 2 lc rgb 'black' axes x1y2 t 'sites'
set size 0.25, 0.347
set rmargin 8
unset tmargin
set origin 0.77,0.61
set y2tics format '%g'
set style histogram rows
set style fill solid 
set boxwidth 0.8
set autoscale x
set autoscale y2
set y2range [0:]
unset ylabel
unset ytics
set y2tics font '/Library/Fonts/Arial Rounded Bold.ttf,8' 
unset grid
unset xlabel
set title 'Expression' offset 0,-1 font '/Library/Fonts/Arial Rounded Bold.ttf,11'
set xrange [] noreverse
set y2label 'total' font '/Library/Fonts/Arial Rounded Bold.ttf,9' offset -1
set style rectangle fc rgb \"black\" 
unset xtics
plot \"plots/$1_exp_tot.txt\" using 2:xticlabels(1) with boxes lt -1 axes x1y2
set origin 0.77,0.33
unset title
set y2label 'full-length' font '/Library/Fonts/Arial Rounded Bold.ttf,9' offset -1
set style fill solid 1.0
set size 0.25, 0.317
plot \"plots/$1_exp_fl.txt\" using 2:xticlabels(1) with boxes lt 3 axes x1y2
set origin 0.77,0.018
set y2label 'optimal' font '/Library/Fonts/Arial Rounded Bold.ttf,9' offset -1
set size 0.25, 0.347
set xtics font '/Library/Fonts/Arial Rounded Bold.ttf,8' offset 0,0.5
plot \"plots/$1_exp_use.txt\" using 2:xticlabels(1) with boxes lt 1 axes x1y2
unset multiplot" > plots/$1_script.txt
cat plots/$1_script.txt | gnuplot
open plots/$1_4nqo.png
#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites',\
#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites'
#  \"$1_edge.txt\" using 1:5 with l lt 4 lw 1 t 'cumulative diff'

#     \"plots/$1.txt\" using 2:29 with l lc rgb 'red' lw 2 t ' DH STD',\
#     \"plots/$1.txt\" using 2:30 with l lc rgb 'pink' lw 2 t ' DH Caff',\
#     \"plots/$1.txt\" using 2:27 with l lc rgb 'blue' lw 2  t '  Ben STD',\
#     \"plots/$1.txt\" using 2:28 with l lc rgb 'cyan' lw 2  t '  Ben YRA1',\

