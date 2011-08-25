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
set output 'figs/density-1D.eps'

set key noauto inside top
#set title 'Density'

# set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'z (nm)'
set ylabel 'density (nm^{-3})'

set style line 1 lt 1 lw 3
set style line 2 lt 1 lc 7 lw 2
set style line 3 lt 3 lw 3


set arrow from 0,0 to 0,0.006 nohead lw 2
set arrow from 4,0 to 4,0.006 nohead lw 2

nl=0.004938863        # liquid density    
nm = 18.8972613       # 1 nm equals this many bohrs

plot [:] [:] \
'figs/cavitysize-04.0.dat' u ($3/nm-2):4 notitle with lines ls 1 , \
'figs/cavitysize-04.0.dat' u ($3/nm-2):(nl) notitle with lines ls 3