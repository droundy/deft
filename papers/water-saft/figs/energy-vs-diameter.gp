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

set terminal postscript eps enhanced dl 3 color "Helvetica" 20
set output 'figs/energy-vs-diameter.eps'

set key noauto outside top

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'Radius (nm)'
set ylabel '{/Symbol m}_{ex}/Area (mN/m)'

set style line 1 lt 1 lc 3 lw 3 pt 7 ps 2
set style line 2 lt 2 lw 3 lc 1 

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin
nl=0.004938863

plot [0:1.0] [0:] \
'figs/single-rod-in-water.dat' u ($1/2):(72.098) notitle with lines lt 0,\
'figs/single-rod-in-water.dat' u ($1/2):($2/(pi*$1*nm)/mNpermeter) notitle with lp ls 1

set size 0.49,0.37        # The second one (inset)
set origin 0.3,0.25
set xlabel
set ylabel
set ytics 2.0
set xtics 0.025
plot [0:0.08] [0:] \
'figs/single-rod-in-water.dat' u ($1/2):($2/(pi*$1*nm)/mNpermeter) notitle with points ls 1,\
'figs/single-rod-in-water.dat' u ($1/2):(kB*298*10*nl*$1/2/mNpermeter) notitle with lines ls 2
