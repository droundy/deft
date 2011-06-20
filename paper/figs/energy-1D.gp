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
set output 'figs/energy-1D.eps'

set key noauto inside top
#set title 'Energy density'

# set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'z (bohr)'
set ylabel 'energy density (Hartree/bohr^{3})'

set style line 1 lt 1 lw 1
set style line 2 lt 1 lc 7 lw 1

plot [:150] [:] \
'figs/cavitysize-060.dat' u 3:5 title 'cavity size = 60 bohr' with lines ls 1 , \
'figs/cavitysize-060.dat' u (45):5 notitle with lines ls 2 , \
'figs/cavitysize-060.dat' u (105):5 notitle with lines ls 2