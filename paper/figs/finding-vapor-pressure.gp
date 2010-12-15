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
set output 'figs/finding-vapor-pressure.eps'

set key noauto

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'density (bohr^{-3})'
set ylabel 'energy density (Hartree/bohr^3)'

nl=4.93889420e-03

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3

plot [:] [:] \
'figs/finding-vapor-pressure.dat' u 1:2 notitle with lines ls 1, \
'figs/finding-vapor-pressure.dat' u 1:3 notitle with lines ls 2

set size 0.38,0.4        # The second one (inset)
set origin 0.58,0.5

#set logscale x          # New command to change x-axis to log-scale

#set xlabel "x (logarithmic scale)" font "Helvetica,16" 0,0.8
unset xlabel
unset ylabel
set xtics 1e-7 font "Helvetica,16"
set ytics 1e-10 font "Helvetica,16"

plot [:3e-7] [:] \
'figs/finding-vapor-pressure.dat' u 1:2 notitle with lines ls 1, \
'figs/finding-vapor-pressure.dat' u 1:3 notitle with lines ls 2





# set title ''

# #set format y "%3.00g"
# set ylabel 'pressure (Hartree/bohr^3)'
# #set ytics 1.0

# #set ticscale 2 1
# set samples 10000

# set xlabel 'temperature (K)'
# #set xtics 50.0

# nl=4.93889420e-03

# set style line 1 lt 1 lw 3
# set style line 2 lt 2 lw 3
# set style line 3 lt 3 lw 3


# plot [:] [:] \
# 'figs/finding-vapor-pressure.dat' u 1:2 notitle with lines ls 1, \
# 'figs/finding-vapor-pressure.dat' u 1:3 notitle with lines ls 2
