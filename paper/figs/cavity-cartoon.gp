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

set terminal postscript eps enhanced color solid lw 3 "Helvetica" 20 size 3,1
set output 'figs/cavity-cartoon.eps'

set key noauto inside top

# set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
#set xlabel 'z (bohr)'
#set ylabel 'pressure (Hartree/bohr^{3})'

unset xtics
unset ytics

set style line 1 lt 3 lw 2

# In cavitysize-xxx.dat columns are 1:x, 2:y, 3:z, 4:density, 5:energy density, 6:entropy, 7:xassoc

plot [45:105] [:0.045] \
'figs/cavitysize-060.dat' u 3:6 notitle with lines ls 1
