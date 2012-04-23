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
set output 'figs/xassoc-single-rod.eps'

set key noauto inside bottom
#set title 'Xassociation'

# set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'Radius (nm)'
set ylabel 'Bonds per Molecule'

set style line 1 lt 1 lw 3 pt 7
set style line 2 pt 6 ps 2 lc 7

set arrow from 0.2,0 to 0.2,4 nohead lt 3 lc 7 lw 1
#set arrow from 4,0 to 4,4 nohead lw 2

nm = 18.8972613       # 1 nm equals this many bohrs

plot [:1] [0:4] \
'figs/single-rod-slice-00.4.dat' u ($2/nm):(4*(1-$7)) notitle with l ls 1 