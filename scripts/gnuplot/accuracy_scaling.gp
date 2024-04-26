
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

set format x '10^{%T}'

set xlabel 'n'
set ylabel 'Relative error'

#set xtics x_0, dx, x_n
#set ytics y_0, dy, y_n

#set cbrange [0:1]

df1 = "../output/analyticsdf_errors.csv"
df2 = "../output/rough_naive_sdf_errors.csv"
df3 = "../output/rough_split_sdf_errors.csv"

set logscale x
set logscale y
set logscale y2

set key top left

set yrange  [1e-16:0.025]
set y2range [1e-16:0.025]

set xrange [900:100000]

#set xtics (1000, 3000, 10000, 30000, 100000)
set ytics (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)

set output "../output/figures/errors.pdf"

set multiplot layout 1,3 margins .1,0.9,0.15,0.9 spacing .05,0.025

  set key font ",8"
  set key maxrows 4 center
  set title "S_1"
  plot df1 using 1:2 with linespoints pointtype 7 pointsize 0.5 lc rgb C1 title "(debiased)", \
       df1 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C2 title "r = 0", \
       df1 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
       df1 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
       df1 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
       df1 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

  unset key
  unset ylabel
  set format y ""
  set title "S_2,  no splitting at 0"
  plot df2 using 1:2 with linespoints pointtype 7 pointsize 0.5 lc rgb C1 title "(debiased)", \
       df2 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C2 title "r = 0", \
       df2 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
       df2 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
       df2 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
       df2 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

  unset key
  set y2tics
  set y2tics (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
  set title "S_2,  split at 0"
  plot df3 using 1:2 with linespoints pointtype 7 pointsize 0.5 lc rgb C1 title "(debiased)", \
       df3 using 1:3 with linespoints pointtype 7 pointsize 0.5 lc rgb C2 title "r = 0", \
       df3 using 1:4 with linespoints pointtype 7 pointsize 0.5 lc rgb C3 title "r = 2", \
       df3 using 1:5 with linespoints pointtype 7 pointsize 0.5 lc rgb C4 title "r = 32", \
       df3 using 1:6 with linespoints pointtype 7 pointsize 0.5 lc rgb C5 title "r = 64", \
       df3 using 1:7 with linespoints pointtype 7 pointsize 0.5 lc rgb C6 title "r = 128"

unset multiplot

unset output

