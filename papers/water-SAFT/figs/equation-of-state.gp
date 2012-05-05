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
set output 'figs/equation-of-state.eps'

set title ''

set key noauto

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0

#set format y "%3.00g"
set ylabel 'Vapor Pressure (atm)'
#set ytics 1.0

#set ticscale 2 1
#set samples 10000

set xlabel 'Temperature (K)'
#set xtics 50.0
set xtics 100
set ytics 25

nl=4.93889420e-03

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 3

atmospheric_pressure = 101325*3.3989316e-14

plot [:] [:100] \
'figs/equation-of-state.dat' u 1:($2/atmospheric_pressure) title 'theory' with lines ls 1, \
'figs/experimental-equation-of-state.dat' u 1:($2/atmospheric_pressure) title 'experiment' with lines ls 2

set size 0.6,0.6        # The second one (inset)
set origin 0.11,0.22
unset xlabel
unset ylabel
set xtics 20 font "Helvetica,16"
#set ytics 1e-9 font "Helvetica,16"
set ytics 0.5 font "Helvetica,16"

plot [273:373] [:] \
'figs/equation-of-state.dat' u 1:($2/atmospheric_pressure) notitle with lines ls 1, \
'figs/experimental-equation-of-state.dat' u 1:($2/atmospheric_pressure) notitle with lines ls 2
