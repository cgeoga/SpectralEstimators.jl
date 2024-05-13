
set terminal cairolatex pdf color size 13cm,10cm
set datafile separator ","

set key off

p1f="../output/p1_fits.csv"
p2f="../output/p2_fits.csv"

set output "fits.tex"

set multiplot layout 2,3 margins .1,0.9,0.15,0.9 spacing .05,0.2

  set yrange [1.5:8.0]
  set xrange [1.5:8.0]
  #set xrange [4.5:5.5]

  set xlabel "$\\hat{\\theta}_1^{\\tiny \\textrm{MLE}}$"
  set ylabel "$\\hat{\\theta}_1$"

  set arrow from 1.5,1.5 to 8.0,8.0 nohead lc rgb "black" lt 2 dt 3
  set arrow from 5.0,1.5 to 5.0,8.0 nohead lc rgb "black" lt 2 dt 3
  set arrow from 1.5,5.0 to 8.0,5.0 nohead lc rgb "black" lt 2 dt 3
  
 
  set title "Whittle" 
  plot p1f using 1:2 title "whittle" pointtype 7 pointsize 0.5 lc rgb "black"

  unset ylabel
  set format y ''
  set title "Debiased Whittle"
  plot p1f using 1:3 title "deb. whittle" pointtype 7 pointsize 0.5 lc rgb "black"


  set title "Our estimator"
  plot p1f using 1:4 title "our estimator" pointtype 7 pointsize 0.5 lc rgb "black"


  unset format y
  set yrange [23.5:40.5]
  set xrange [23.5:40.5]
  #set xrange [34.5:35.5]

  set ytics (24,28,32,36,40)

  unset title

  set xtics (24,28,32,36,40)

  set xlabel "$\\hat{\\theta}_2^{\\tiny \\textrm{MLE}}$"
  set ylabel "$\\hat{\\theta}_2$"

  set arrow from 23.5,23.5 to 40.5,40.5 nohead lc rgb "black" lt 2 dt 3
  set arrow from 35.0,23.5 to 35.0,40.5 nohead lc rgb "black" lt 2 dt 3
  set arrow from 23.5,35.0 to 40.5,35.0 nohead lc rgb "black" lt 2 dt 3

  plot p2f using 1:2 title "whittle" pointtype 7 pointsize 0.5 lc rgb "black"

  unset ylabel
  set format y ''

  plot p2f using 1:3 title "deb. whittle" pointtype 7 pointsize 0.5 lc rgb "black"
  plot p2f using 1:4 title "our estimator" pointtype 7 pointsize 0.5 lc rgb "black"

unset multiplot

unset output

