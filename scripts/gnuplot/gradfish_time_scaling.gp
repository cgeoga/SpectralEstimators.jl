
set terminal pdf color size 16cm,6cm enhanced font 'Verdana,10'
#set terminal epslatex color blacktext size 16cm,8cm
set datafile separator ","
set autoscale xfix
set autoscale yfix

#set cbtics format '\tiny %g'
#set cblabel '\tiny label'

#set format x '\tiny %g'
#set format y '\tiny %g'

load "colors.gp"

set xlabel 'n'
set ylabel 'Time (s)'

set format x '10^{%T}'

#set xtics x_0, dx, x_n
#set ytics y_0, dy, y_n

#set cbrange [0:1]

df1 = "../output/rough_sdf_grad_times.csv"
df2 = "../output/rough_sdf_exact_fish_times.csv"
df3 = "../output/rough_sdf_stoch_fish_times.csv"

set logscale x
set logscale y
set logscale y2

set key top left

set xrange [900:50000]
set yrange [0.025:50.0]
set y2range [0.025:50.0]

set ytics (0.1, 1.0, 10.0, 40.0)

set output "../output/figures/gradfish_runtime.pdf"

set multiplot layout 1,3 margins .085,0.915,0.2,0.9 spacing .025,0.025

unset key
set title "gradient"
plot df1 using 1:2 with lines dashtype 2 linewidth 4 lc rgb C7 title "O(n log n)", \
     df1 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
     df1 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C8 title "r = 4", \
     df1 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
     df1 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
     df1 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

unset ylabel
set format y ''
set format y2 ''
set y2tics (0.1, 1.0, 10.0, 40.0)
set title "exact expected Fisher"
plot df2 using 1:2 with lines dashtype 2 linewidth 4 lc rgb C7 title "O(n log n)", \
     df2 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
     df2 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C8 title "r = 4", \
     df2 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
     df2 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
     df2 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

unset format y2
set y2tics (0.1, 1.0, 10.0, 40.0)
set title "stochastic expected Fisher"
set key font ",8"
set key maxrows 4 top left
plot df3 using 1:2 with lines dashtype 2 linewidth 4 lc rgb C7 title "O(n log n)", \
     df3 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
     df3 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C7 title "r = 4", \
     df3 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
     df3 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
     df3 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

unset multiplot

unset output


