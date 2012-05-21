#!/usr/bin/gnuplot
#
# guide for line and point styles:
#
#  0  ..............  .                    broken line
#  1  --------------  +                    red
#  2  -- -- -- -- --  x                    green
#  3  -  -  -  -  -   *                    blue
#  4  ..............  empty square         magenta
#  5  __.__.__.__.__  full  square         cyan
#  6  _ . _ . _ . _   empty circle         yellow
#  7  - -  - -  - -   full  circle         black
#  8  - - -  - - -    empty up triangle    brown
#  9  - - - -  - - -  full  up triangle    grey
# 10 (1)              empty down triangle
# 11 (2)              full  down triangle
# 12 (3)              empty diamond
# 13 (4)              full  diamond
# 14 (5)              empty pentagon
# 15 (6)              full  pentagon
# 16-31               watches

set terminal postscript eps enhanced color "Helvetica" 20
set output 'figs/density-single-rod.eps'

reset
unset arrow
set view map

set key inside bottom
#set title 'One hydrophobic rod - Density'

# set multiplot
#set pm3d map
#set palette color positive
#set ticslevel 0
#set samples 50; set isosamples 50
#set palette rgbformulae 22,13,-31

#set size 1,1          # The first plot (host plot)
#set origin 0,0
#set xlabel 'y (nm)'
#set ylabel 'z (nm)'

set style line 1 lt 1 lc 1 lw 3
set style line 2 lt 1 lc 3 lw 3
set style line 3 lt 1 lc 7 lw 3
set style line 4 lt 3 lc 3 lw 3

#splot [:] [:] [:] \
#'figs/single-rod-01.0.dat' u ($2/nm):($3/nm):4 notitle with pm3d 

set xlabel 'Radius (nm)'
set ylabel 'Density (g/mL)'

set label "r=0.3nm" at 0.62,1.74 rotate by 0 font 'Helvetica,20' 
set arrow from 0.6,1.72 to 0.46,1.6 lw 2
set label "r=0.5nm" at 0.75,1.38 rotate by 0 font 'Helvetica,20' textcolor lt 3
set arrow from 0.73,1.36 to 0.55,1.22 lw 2 lc rgb "blue"
set label "r=1.0nm" at 1.25,0.75 rotate by 0 font 'Helvetica,20' textcolor lt 1
set arrow from 1.23,0.74 to 1.04,0.6 lw 2 lc rgb "red"

nl=0.004938863
nm = 18.8972613     # 1 nm equals this many bohrs
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

plot [:1.50] [:] \
'figs/single-rod-slice-02.0.dat' u ($2/nm):($4/gpermL) notitle with lines ls 1 , \
'figs/single-rod-slice-01.6.dat' u ($2/nm):($4/gpermL) notitle with lines ls 2 , \
'figs/single-rod-slice-01.0.dat' u ($2/nm):($4/gpermL) notitle with lines ls 2 , \
'figs/single-rod-slice-00.6.dat' u ($2/nm):($4/gpermL) notitle with lines ls 3 , \
'figs/single-rod-slice-00.3.dat' u ($2/nm):($4/gpermL) notitle with lines ls 3 , \
'figs/single-rod-slice-00.1.dat' u ($2/nm):($4/gpermL) notitle with lines ls 1 , \
'figs/single-rod-slice-00.6.dat' u ($2/nm):(nl/gpermL) notitle with lines ls 4 , \
'figs/single-rods-calculated-density.dat' u ($1/nm):($2/gpermL) notitle with points ls 1
