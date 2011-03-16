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
set output 'figs/pressure-with-isotherms.eps'

set key noauto

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'density (bohr^{-3})'
set ylabel 'pressure ( )'

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3

kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

plot [:] [-2e-6:2e-6] \
'figs/pressure-with-isotherms.dat' u 1:3 title 'pressure at T = 293K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:5 title 'pressure at T = 343K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:7 title 'pressure at T = 393K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:9 title 'pressure at T = 443K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:11 title 'pressure at T = 493K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:13 title 'pressure at T = 543K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:15 title 'pressure at T = 593K' with lines , \
'figs/pressure-with-isotherms.dat' u 1:17 title 'pressure at T = 643K' with lines , \
'figs/equation-of-state.dat' u 3:2 title 'eos pressure' with lines ls 1 , \
'figs/equation-of-state.dat' u 4:2 notitle with lines ls 1 , \
'figs/experimental-equation-of-state.dat' u 3:2 title 'expt pressure' with lines ls 2 , \
'figs/experimental-equation-of-state.dat' u 4:2 notitle with lines ls 2 , \
'figs/pressure-with-isotherms.dat' u 1:(0) notitle with lines 