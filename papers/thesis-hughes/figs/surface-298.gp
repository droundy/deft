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
set output 'figs/surface-298.eps'

set key noauto

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set ylabel 'Density (g/mL)'
set xlabel 'Position (nm)'

nl=0.004938863
nm = 18.8972613       # 1 nm equals this many bohrs
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

set style line 1 lt 1 lw 3
set style line 2 lt 3 lc 3 lw 3

plot [:4] [:] \
'figs/surface-298.dat' u ($1/nm):($2/gpermL) notitle with lines ls 1, \
'figs/surface-298.dat' u ($1/nm):(0.996782051315) notitle with lines ls 2
