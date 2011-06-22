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
set output 'figs/rods-energy-vs-distance.eps'

set key noauto inside top
#set title 'Energy of 2 rods vs distance between them'

# set multiplot

set encoding iso_8859_1

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'd (nm)'
set ylabel 'Energy per length (kJ/mol {\305})'

set style line 1 lt 1 pt 7 ps 2 lw 2

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
eV = 0.036749326 # Hartree
kJpermol = 1.04e-2*eV # in atomic units
angstrom = 0.1*nm

plot [:] [:] \
'figs/rods-in-water.dat' u 1:($2/(kJpermol/angstrom) - 59.06) notitle with p ls 1

#'figs/rods-in-water.dat' u 1:($2/(2*pi*1*nm)/mNpermeter/2) notitle with lines ls 1 
