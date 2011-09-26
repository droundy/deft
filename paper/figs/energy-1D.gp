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
set xlabel 'z (nm)'
set ylabel 'Energy Density (kJ/mol nm^{3})'

set style line 1 lt 1 lc 3 lw 2
set style line 2 lt 1 lc 7 lw 1
set style line 3 lt 1 lc 2 lw 2
set style line 4 lt 1 lc 1 lw 2


set arrow from 0,-4e-5 to 0,1.1e-5 nohead lw 2
set arrow from 4,-4e-5 to 4,1.1e-5 nohead lw 2

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
eV = 0.036749326 # Hartree
kJpermol = 1.04e-2*eV # in atomic units
angstrom = 0.1*nm
kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

plot [:] [:] \
'figs/cavitysize-04.0.dat' u ($3/nm-2):($5/kJpermol*nm**3) title 'F' with lines ls 1 , \
'figs/cavitysize-04.0.dat' u ($3/nm-2):(($5-kB*298*$6)/kJpermol*nm**3) title 'U' with lines ls 4 , \
'figs/cavitysize-04.0.dat' u ($3/nm-2):((-kB*298*$6)/kJpermol*nm**3) title '-TS' with lines ls 3 
