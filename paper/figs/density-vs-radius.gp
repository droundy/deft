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

set terminal postscript eps enhanced color solid "Helvetica" 20
set output 'figs/density-vs-radius.eps'

reset
unset arrow
set view map

set key inside top

set xlabel 'radius (nm)'
set ylabel 'density (bohr^{-3})'
set ytics 0.0005
set style line 1 lt 7 lw 1

set style line 2 pt 6 lc 7
set style line 3 lt 1 lw 2
set style line 4 pt 1 lc 3

set pointsize 1.5

nl=0.004938863
nm = 18.8972613
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

plot [0.5:1.5] [0.004:] \
'figs/single-rod-res0.05-slice-01.0.dat' u ($3/nm):4 title 'high resolution' with lines ls 3 , \
'figs/single-rod-slice-01.0.dat' u ($2/nm):4 title 'medium resolution' with points ls 2 , \
'figs/single-rod-res0.5-slice-01.0.dat' u ($2/nm):4 title 'low resolution' with points ls 4
