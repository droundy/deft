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

set terminal postscript eps enhanced color "Helvetica" 26
set output 'figs/surface-tension.eps'

set title ''

#set format y "%3.00g"
set ylabel 'Surface Tension (mN/m)'
#set ytics 1.0

#set ticscale 2 1
set samples 10000

set xlabel 'Temperature (^{o}C)'
#set xtics 10.0

mNpermeter = 6.4230498e-07 # in atomic units

nl=4.93889420e-03

set style line 1 lt 1 lw 3
set style line 2 lt 3 lw 3
set style line 3 lt 3 lw 3

set multiplot

set size 1,1          # The first plot (host plot)
plot [0:] [:] \
'figs/surface-tension.dat' u ($1-273.15):($2/mNpermeter) title 'theory' with lines ls 1, \
'figs/experimental-surface-tension.csv' u ($1-273.15):($2/mNpermeter) title 'experiment' with lines ls 2

set size 0.45,0.35        # The second one (inset)
set origin 0.12,0.2
set xlabel
set ylabel
set ytics 2.0
plot [0:50] [68:76] \
'figs/surface-tension.dat' u ($1-273.15):($2/mNpermeter) notitle with lines ls 1, \
'figs/experimental-surface-tension.csv' u ($1-273.15):($2/mNpermeter) notitle with lines ls 2
