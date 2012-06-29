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
set output 'figs/rods-energy-vs-distance.eps'

set key at 1.2,-0.18 # inside bottom
#set title 'Energy of 2 rods vs distance between them'

#set multiplot

#set encoding iso_8859_1

#set size 1,1          # The first plot (host plot)
#set origin 0,0
set xlabel 'd (nm)'
set ylabel 'Free energy per length (kJ/mol nm)'
set mxtics 2

set style line 1 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#eecc11"
set style line 2 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#aacc11"
set style line 3 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#11cc33"
set style line 4 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#008888"
set style line 5 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#002266"

set style line 11 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#eecc11"
set style line 12 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#aacc11"
set style line 13 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#11cc33"
set style line 14 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#008888"
set style line 15 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#002266"

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
eV = 0.036749326 # Hartree
kJpermol = 1.04e-2*eV # in atomic units
angstrom = 0.1*nm
kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

#'figs/rods-in-water.dat' u 1:($2/(2*pi*1*nm)/mNpermeter/2) notitle with lines ls 1 
#'figs/rods-in-water.dat' u 1:($2/(kJpermol/angstrom) - 59.06) notitle with lp ls 1

set label "Before transition" at 0.05,0.03 rotate by 0 font 'Helvetica,20'
#set label "(FIG. 14 top)" at 0.05,0.017 rotate by 0 font 'Helvetica,20' 
set arrow from 0.16,0.022 to 0.195,-0.075 lw 2
set label "After transition" at 0.75,0.03 rotate by 0 font 'Helvetica,20' 
#set label "(FIG. 14 bottom)" at 0.215,0.017 rotate by 0 font 'Helvetica,20'
set arrow from 0.74,0.03 to 0.61,0.0075 lw 2

gamma = 72 # mN/m
f(d,r) = d>(pi-2)*r/2 ? 0 : (d - (pi-2)*r/2)*2*gamma*mNpermeter/kJpermol

plot [:1.3] [:0.05] \
'figs/rods-in-water-00.6nm.dat' u 1:($2/kJpermol/(nm)-0.411) title 'r=0.3 nm' with l ls 1 , \
'figs/rods-in-water-00.6nm.dat' u 1:(f($1,0.6)) notitle with l ls 11 , \
'figs/rods-in-water-01.0nm.dat' u 1:($2/kJpermol/(nm)-0.776) title 'r=0.5 nm' with l ls 2 , \
'figs/rods-in-water-01.0nm.dat' u 1:(f($1,1.0)) notitle with l ls 12 , \
'figs/rods-in-water-01.4nm.dat' u 1:($2/kJpermol/(nm)-1.127) title 'r=0.7 nm' with l ls 3 , \
'figs/rods-in-water-01.4nm.dat' u 1:(f($1,1.4)) notitle with l ls 13 , \
'figs/rods-in-water-01.8nm.dat' u 1:($2/kJpermol/(nm)-1.477) title 'r=0.9 nm' with l ls 4 , \
'figs/rods-in-water-01.8nm.dat' u 1:(f($1,1.8)) notitle with l ls 14 , \
'figs/rods-in-water-02.0nm.dat' u 1:($2/kJpermol/(nm)-1.652) title 'r=1.0 nm' with l ls 5, \
'figs/rods-in-water-02.0nm.dat' u 1:(f($1,2.0)) notitle with l ls 15
#'figs/rods-in-water-02.4nm.dat' u 1:($2/kJpermol/(nm)-1.939) title 'r=1.2 nm' with lp ls 6
