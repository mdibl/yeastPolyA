grep -E "$1$"  saccer3_simple_features.txt | cut -f 3,4 | perl -ane "print \"@F[0]\t-0.03\n@F[1]\t-0.03\n\";" > plots/$1_coord.txt 
grep -E "$1[[:space:]]" dh_all_Arem_1000_05_plots.txt | grep -v totals | grep -v bounds | grep -v DH > plots/$1.txt
SEN=`grep -E "$1$"  saccer3_simple_features.txt | cut -f 5`
if [ $SEN == "-" ]; then
   echo plotting antisense
   REV="reverse"   
fi

#grep tota plots/$1.txt | cut -f 3-4 | perl -ane "print \"STD\t\"; print (1*@F[1]/476); print \"\n4NQO\t\"; print (1*@F[0]/1363); print \"\n\"; " > plots/$1_exp.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 2-3 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_tot.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 4-5 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_fl.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 6-7 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_use.txt
#./efTest plots/$1.txt 1 3 12 0.05 100 > plots/$1_edge.txt
#echo "set term aqua font 'Arial,16' size 900 600 title '$1' 
echo "set term png size 800,640 font \"/Library/Fonts/Arial Narrow.ttf,10\" 
set output 'plots/$1_all.png'
set multiplot
set key inside top center horizontal box
set size 0.95, 0.65
set origin 0.025, 0.35
set grid lt 0
set title '$1: `grep -E $1[[:space:]] alias_commas.txt | cut -f 2`' offset 1.0 font '/Library/Fonts/Arial Narrow.ttf,16'
set ylabel 'cumulative poly(A) probability' font '/Library/Fonts/Arial Narrow.ttf,14' 
set xtics font '/Library/Fonts/Arial Narrow.ttf,11' 
set ytics font '/Library/Fonts/Arial Narrow.ttf,11' 
set xrange [] $REV writeback
set xtics format ''
set yrange [-0.1:1.15]
set bmargin 1
plot \"plots/$1.txt\" using 2:3 with l lc rgb 'black' lw 2 t 'WT',\
     \"plots/$1.txt\" using 2:5 with l lc rgb 'blue'  lw 2 t 'pcf11-13',\
     \"plots/$1_coord.txt\" using 1:2 with l lt 9 lw 8 t 'CDS'
unset title
set xtics format \"%.0f\" font '/Library/Fonts/Arial Narrow.ttf,12'
set xlabel 'position (nt)' font '/Library/Fonts/Arial Narrow.ttf,14' 
set ylabel 'difference'
set yrange [-0.7:0.7]
set xrange restore
set size 0.95,0.35
set origin 0.025,0.0
unset bmargin
set tmargin 0
plot 'plots/$1.txt' u 2:(\$5-\$3) w l lc rgb 'black' lw 2 t 'difference'
unset multiplot" > plots/$1_script.txt
cat plots/$1_script.txt | gnuplot
open plots/$1_all.png
#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites',\
#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites'
#  \"$1_edge.txt\" using 1:5 with l lt 4 lw 1 t 'cumulative diff'
#
