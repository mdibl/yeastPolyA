grep -E "$1$"  saccer3_simple_features.txt | cut -f 3,4 | perl -ane "print \"@F[0]\t-0.03\n@F[1]\t-0.03\n\";" > plots/$1_coord.txt 
grep -E "$1[[:space:]]" dh_all_Arem_1000_05_plots.txt | grep -v totals | grep -v bounds | grep -v DH > plots/$1.txt
SEN=`grep -E "$1$"  saccer3_simple_features.txt | cut -f 5`
ARR1="head"
ARR2="tail"
if [ $SEN == "-" ]; then
   echo plotting antisense
   REV="reverse"
   ARR1="tail"
   ARR2="head"
fi
#sh makeExpressionGrid.sh $1
sh makeExpressionGrid_tx.sh $1
#grep tota plots/$1.txt | cut -f 3-4 | perl -ane "print \"STD\t\"; print (1*@F[1]/476); print \"\n4NQO\t\"; print (1*@F[0]/1363); print \"\n\"; " > plots/$1_exp.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 2-3 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_tot.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 4-5 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_fl.txt
#grep -E "$1[[:space:]]" 4nqo_count_summary.txt | cut -f 6-7 | perl -ane "print \"STD\t@F[1]\n4NQO\t@F[0]\n\"; " > plots/$1_exp_use.txt
#./efTest plots/$1.txt 1 3 12 0.05 100 > plots/$1_edge.txt
#echo "set term aqua font 'Arial,16' size 900,800 title '$1' 
echo "set term png size 1200,800 font \"/Library/Fonts/Arial Narrow.ttf,10\" 
set output 'plots/$1_bothWC.png'
set multiplot
#set key outside top right vertical Right  width 0 spacing 2 box font '/Library/Fonts/Arial Narrow.ttf,9' height 1.25
unset key
set size 0.7, 0.475
set origin 0.025, 0.5
set grid lt 0
set style arrow 3 head filled size screen 0.014,15,45 lc rgb 'blue' lw 4
set arrow from `$ARR1 -n 1 plots/$1_coord.txt | cut -f 1`,`$ARR1 -n 1 plots/$1_coord.txt | cut -f 2` to `$ARR2 -n 1 plots/$1_coord.txt | cut -f 1`,`$ARR2 -n 1 plots/$1_coord.txt | cut -f 2` as 3
set title '$1: `grep -E $1[[:space:]] alias_commas.txt | cut -f 2`' offset 20,-1 font '/Library/Fonts/Arial Narrow.ttf,16'
set ylabel 'cumulative poly(A) probability' font '/Library/Fonts/Arial Narrow.ttf,14' 
set xtics font '/Library/Fonts/Arial Narrow.ttf,11' 
set ytics font '/Library/Fonts/Arial Narrow.ttf,11' 
set xrange [] $REV writeback
set xtics format ''
set yrange [-0.1:1.15]
set bmargin 1
set tmargin 2
set rmargin 19
set label 1 'PCF11-dQN' at graph 0.02,0.93 font '/Library/Fonts/Arial Narrow.ttf,14'
plot \"plots/$1.txt\" using 2:7 with l lc rgb 'black' lw 2 t 'WT',\
     \"plots/$1.txt\" using 2:8 with l lc rgb 'gray'  lw 2 t 'caff + WT',\
     \"plots/$1.txt\" using 2:9 with l lc rgb 'red'  lw 2 t 'dQN',\
     \"plots/$1.txt\" using 2:10 with l lc rgb 'pink'  lw 2 t 'caff + dQN'
unset title
unset bmargin
set tmargin 1
set xtics format \"%.0f\" font '/Library/Fonts/Arial Narrow.ttf,12'
set xrange restore
set origin 0.025,0.025
set label 1 'PCF11-13' at graph 0.02,0.93 font '/Library/Fonts/Arial Narrow.ttf,14'
plot \"plots/$1.txt\" using 2:3 with l lc rgb 'black' lw 2 t 'WT',\
     \"plots/$1.txt\" using 2:4 with l lc rgb 'gray'  lw 2 t 'caff + WT',\
     \"plots/$1.txt\" using 2:5 with l lc rgb 'red'  lw 2 t '11-13',\
     \"plots/$1.txt\" using 2:6 with l lc rgb 'pink'  lw 2 t 'caff + 11-13'
set origin 0.63,0.48
#set boxwidth 0.9
unset ylabel
set bmargin 6
unset tmargin
unset rmargin
set xrange [] noreverse
set autoscale y
set yrange [0:]
set autoscale x
set key outside center bmargin horizontal Right box height 1
set size 0.35,0.44
set style fill solid 1.00 border lt -1
set style histogram clustered gap 1 title textcolor lt -1
set datafile missing '-'
set style data histograms
set xtics border in scale 0,0 nomirror rotate by 0  autojustify  font '/Library/Fonts/Arial Narrow.ttf,10'
unset label 1
set title \"Expression\" offset 0,-1  font '/Library/Fonts/Arial Narrow.ttf,14'
set palette grey
plot \"plots/$1_pcf11-dQN_exp_tx.txt\" using 2:xtic(1) ti col fc rgb 'black', '' u 3 ti col fc rgb 'gray', '' u 4 ti col fc rgb 'red', '' u 5 ti col fc rgb 'pink'
set origin 0.63,0.015
unset key
plot \"plots/$1_pcf11-13_exp_tx.txt\" using 2:xtic(1) ti col fc rgb 'black', '' u 3 ti col fc rgb 'gray', '' u 4 ti col fc rgb 'red', '' u 5 ti col fc rgb 'pink'
unset multiplot" > plots/$1_script.txt
cat plots/$1_script.txt | gnuplot
open plots/$1_bothWC.png

#set key inside right top vertical Right noreverse noenhanced autotitle nobox

#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites',\
#  \"$1_edge.txt\" using 1:4 with l lt 3 lw 2 t 'APA sites'
#  \"$1_edge.txt\" using 1:5 with l lt 4 lw 1 t 'cumulative diff'
#
#plot 'plots/$1.txt' u 2:(\$4-\$3) w l lc rgb 'black' lw 2 t 'WT caff', '' u 2:(\$6-\$5) with l lc rgb 'red' lw 2 t 'pcf11-13 caff'
#set xlabel 'position (nt)' font '/Library/Fonts/Arial Narrow.ttf,14' 
#set yrange [ 0.00000 : 300000. ] noreverse nowriteback
#set ylabel 'difference'
#set yrange [-0.7:0.7]
#set size 0.95,0.35
