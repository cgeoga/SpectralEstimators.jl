
set terminal cairolatex pdf color size 13cm,10cm
set datafile separator ","
set autoscale xfix
set autoscale yfix

load "colors.gp"
load "pointtypes.gp"

#set cbtics format '\tiny %g'
#set cblabel '\tiny label'

set ylabel 'Relative error'

#set xtics x_0, dx, x_n
#set ytics y_0, dy, y_n

#set cbrange [0:1]

set format y  '$10^{%T}$'
set format y2 '$10^{%T}$'

df1 = "../output/analyticsdf_grad_errors.csv"
df2 = "../output/rough_split_sdf_grad_errors.csv"

df3  = "../output/analyticsdf_exact_fish_errors.csv"
df3s = "../output/analyticsdf_stoch_fish_errors.csv"
df4  = "../output/rough_split_sdf_exact_fish_errors.csv"
df4s = "../output/rough_split_sdf_stoch_fish_errors.csv"

set logscale x
set logscale y
set logscale y2

set key center left

set yrange  [1e-16:0.1]
set y2range [1e-16:0.025]

set xrange [900:50000]

set ytics (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)

set output "gradfish_errors.tex"

set multiplot layout 2,2 margins .1,0.9,0.125,0.9 spacing .035,0.15

  unset xlabel 

  unset key

  unset xtics
  set xtics
  unset y2tics
  set format x ''
  set yrange   [1e-16:0.1]
  set y2range  [1e-16:0.1]
  set ytics
  set ytics (1e-15, 1e-12, 1e-9, 1e-6, 1e-3)
  set title "$S_1$ gradient"
  plot df1 using 1:2 with linespoints pointtype P3 pointsize 0.5 lc rgb C3 title "$r = 2$", \
       df1 using 1:3 with linespoints pointtype P8 pointsize 0.5 lc rgb C8 title "$r = 4$", \
       df1 using 1:4 with linespoints pointtype P4 pointsize 0.5 lc rgb C4 title "$r = 32$", \
       df1 using 1:5 with linespoints pointtype P5 pointsize 0.5 lc rgb C5 title "$r = 64$", \
       df1 using 1:6 with linespoints pointtype P6 pointsize 0.5 lc rgb C6 title "$r = 128$"

  unset ylabel
  set format y ''
  set y2tics
  set y2tics (1e-15, 1e-12, 1e-9, 1e-6, 1e-3)
  set title "$S_2$ (split at 0) gradient"
  plot df2 using 1:2 with linespoints pointtype P3 pointsize 0.5 lc rgb C3 title "$r = 2$", \
       df2 using 1:3 with linespoints pointtype P8 pointsize 0.5 lc rgb C8 title "$r = 4$", \
       df2 using 1:4 with linespoints pointtype P4 pointsize 0.5 lc rgb C4 title "$r = 32$", \
       df2 using 1:5 with linespoints pointtype P5 pointsize 0.5 lc rgb C5 title "$r = 64$", \
       df2 using 1:6 with linespoints pointtype P6 pointsize 0.5 lc rgb C6 title "$r = 128$"


  set ylabel "Relative error"
  set format x '$10^{%T}$'
  set xlabel "n"
 
  set key font ",6" 
  set key width 9
  set key spacing 1.3
  set key maxrows 4 center at 4000,8e-10
  unset y2tics
  set format y  '$10^{%T}$'
  set ytics (1e-15, 1e-12, 1e-9, 1e-6, 1e-3)
  set title "$S_1$ expected Fisher"
  plot df3 using 1:2  with linespoints pointtype P3  pointsize 0.5 lc rgb C3  title "$r = 2$", \
       df3 using 1:3  with linespoints pointtype P8  pointsize 0.5 lc rgb C8  title "$r = 4$", \
       df3 using 1:4  with linespoints pointtype P4  pointsize 0.5 lc rgb C4  title "$r = 32$", \
       df3 using 1:5  with linespoints pointtype P5  pointsize 0.5 lc rgb C5  title "$r = 64$", \
       df3 using 1:6  with linespoints pointtype P6  pointsize 0.5 lc rgb C6  title "$r = 128$", \
       df3s using 1:2 with linespoints pointtype P9  pointsize 0.5 lc rgb C9  title "$r = 2$ (stoch)", \
       df3s using 1:6 with linespoints pointtype P10 pointsize 0.5 lc rgb C10 title "$r = 128$ (stoch)"


  unset key
  unset ylabel
  set format y ''
  set y2tics (1e-15, 1e-12, 1e-9, 1e-6, 1e-3)
  set title "$S_2$ (split at 0) expected Fisher"
  plot df4 using 1:2  with linespoints pointtype P3  pointsize 0.5 lc rgb C3  title "$r = 2$", \
       df4 using 1:3  with linespoints pointtype P8  pointsize 0.5 lc rgb C8  title "$r = 4$", \
       df4 using 1:4  with linespoints pointtype P4  pointsize 0.5 lc rgb C4  title "$r = 32$", \
       df4 using 1:5  with linespoints pointtype P5  pointsize 0.5 lc rgb C5  title "$r = 64$", \
       df4 using 1:6  with linespoints pointtype P6  pointsize 0.5 lc rgb C6  title "$r = 128$", \
       df4s using 1:2 with linespoints pointtype P9  pointsize 0.5 lc rgb C9  title "$r = 2$ (stoch)", \
       df4s using 1:6 with linespoints pointtype P10 pointsize 0.5 lc rgb C10 title "$r = 128$ (stoch)"

unset multiplot

unset output

