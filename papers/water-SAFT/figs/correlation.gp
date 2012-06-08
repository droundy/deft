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

set terminal postscript landscape enhanced color "Helvetica" 26
set output 'figs/correlation.eps'

set title ''

#set format y "%3.00g"
set ylabel 'n/n_l'
set ytics 1.0

#set ticscale 2 1
set samples 10000

set xlabel 'distance (Angstrom)'
set xtics 1.0

nl=4.93889420e-03

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 3


plot [:5] [:4] \
'figs/gOO.exp' u ($1*0.529):($2+0.015444) title 'g_{OO}' with lines ls 1, \
'figs/gHH.exp' u ($1*0.529):($2+0.00846079) title 'g_{HH}' with lines ls 2, \
'figs/gOH.exp' u ($1*0.529):($2+0.00406799) title 'g_{OH}' with lines ls 3
