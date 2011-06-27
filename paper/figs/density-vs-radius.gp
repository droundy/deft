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
set output 'figs/density-vs-radius.eps'

reset
unset arrow
set view map

set key inside bottom
set title 'One hydrophobic rod - Density'

# set multiplot
#set pm3d map
#set palette color positive
#set ticslevel 0
#set samples 50; set isosamples 50
#set palette rgbformulae 22,13,-31

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'z (nm)'
set ylabel 'density (bohr^{-3})'

set style line 1 lt 7 lw 1
#set style line 2 lt 1 lc 7 lw 1
set style line 2 pt 6
set style line 3 pt 0 lc 4 lw 2

nm = 18.8972613

plot [:2.5] [:] \
'figs/single-rod-slice-01.0.dat' u ($3/nm):4 notitle with points ls 2 , \
'figs/single-rod-res0.05-slice-01.0.dat' u ($3/nm):4 notitle with points ls 3
