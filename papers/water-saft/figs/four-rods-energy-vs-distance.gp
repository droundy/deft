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
set output 'figs/four-rods-energy-vs-distance.eps'

set key at 3.25, -1.1
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
set style line 5 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#002299"
set style line 6 lt 1 pt 7 ps 1.5 lw 5 lc rgb "#002255"

set style line 11 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#eecc11"
set style line 12 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#aacc11"
set style line 13 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#11cc33"
set style line 14 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#008888"
set style line 15 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#002299"
set style line 16 lt 2 pt 7 ps 1.5 lw 5 lc rgb "#002255"

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
eV = 0.036749326 # Hartree
kJpermol = 1.04e-2*eV # in atomic units
angstrom = 0.1*nm
kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

#'figs/rods-in-water.dat' u 1:($2/(2*pi*1*nm)/mNpermeter/2) notitle with lines ls 1 
#'figs/rods-in-water.dat' u 1:($2/(kJpermol/angstrom) - 59.06) notitle with lp ls 1

#set label "Before transition" at 0.05,0.03 rotate by 0 font 'Helvetica,20'
#set label "(FIG. 14 top)" at 0.05,0.017 rotate by 0 font 'Helvetica,20' 
#set arrow from 0.16,0.022 to 0.195,-0.075 lw 2
#set label "After transition" at 0.75,0.03 rotate by 0 font 'Helvetica,20' 
#set label "(FIG. 14 bottom)" at 0.215,0.017 rotate by 0 font 'Helvetica,20'
#set arrow from 0.74,0.03 to 0.61,0.0075 lw 2

gamma = 72.098 # mN/m
f(d,r) = d>(3*pi-4)*r/2 ? 0 : (-6*pi*r+8*r+4*d)*gamma*mNpermeter/kJpermol

plot [:3.5] [:0.1] \
'figs/four-rods-in-water-00.4nm.dat' u 1:($2/kJpermol/(nm)-0.286) title 'r=0.2 nm' with l ls 1, \
'figs/four-rods-in-water-00.4nm.dat' u 1:(f($1,0.2)) notitle with l ls 11 , \
'figs/four-rods-in-water-00.8nm.dat' u 1:($2/kJpermol/(nm)-1.128) title 'r=0.4 nm' with l ls 2, \
'figs/four-rods-in-water-00.8nm.dat' u 1:(f($1,0.4)) notitle with l ls 12 , \
'figs/four-rods-in-water-01.2nm.dat' u 1:($2/kJpermol/(nm)-1.863) title 'r=0.6 nm' with l ls 3, \
'figs/four-rods-in-water-01.2nm.dat' u 1:(f($1,0.6)) notitle with l ls 13 , \
'figs/four-rods-in-water-01.6nm.dat' u 1:($2/kJpermol/(nm)-2.573) title 'r=0.8 nm' with l ls 4, \
'figs/four-rods-in-water-01.6nm.dat' u 1:(f($1,0.8)) notitle with l ls 14 , \
'figs/four-rods-in-water-02.0nm.dat' u 1:($2/kJpermol/(nm)-3.279) title 'r=1.0 nm' with l ls 5, \
'figs/four-rods-in-water-02.0nm.dat' u 1:(f($1,1.0)) notitle with l ls 15, \
'figs/four-rods-in-water-02.4nm.dat' u 1:($2/kJpermol/(nm)-3.988) title 'r=1.2 nm' with l ls 6, \
'figs/four-rods-in-water-02.4nm.dat' u 1:(f($1,1.2)) notitle with l ls 16
