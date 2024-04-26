
set terminal pdf color size 16cm,6cm enhanced font 'Verdana,10'
#set terminal epslatex color blacktext size 16cm,8cm
#set terminal pngcairo size 1200,800 enhanced font 'Verdana,10'
#set key off
set datafile separator ","
set autoscale xfix
set autoscale yfix

load "colors.gp"

#set cbtics format '\tiny %g'
#set cblabel '\tiny label'

#set format x '\tiny %g'
#set format y '\tiny %g'

set xlabel 'n'
set ylabel 'Time (s)'

set format x '10^{%T}'

#set xtics x_0, dx, x_n
#set ytics y_0, dy, y_n

#set cbrange [0:1]

df1 = "../output/rough_naive_sdf_times.csv"

set logscale x
set logscale y

set key font ",9"
set key top left maxrows 4

set xrange [900:100000]


set output "../output/figures/runtime.pdf"

plot df1 using 1:2 with lines linewidth 4 dashtype 2 lc rgb C7 title "O(n log n)", \
     df1 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C1 title "(debiased)", \
     df1 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C2 title "r = 0", \
     df1 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
     df1 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
     df1 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
     df1 using 1:8 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

unset output

