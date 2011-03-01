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
set output 'figs/near-critical-point.eps'

set key noauto

set xlabel 'density (bohr^{-3})'
set ylabel 'energy density (Hartree/bohr^3)'

nl=4.93889420e-03

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3

plot [0.0005:0.003] [:] \
'figs/finding-vapor-pressure.dat' u 1:($4-$5) notitle with lines ls 1
