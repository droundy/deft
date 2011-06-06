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

set key noauto outside bottom
set key font "Helvetica,14"
set key horizontal

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'density (bohr^{-3})'
set ylabel 'pressure (Hartree/bohr^{3})'

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3
set style line 3 pt 6 ps 1 lc 3
set style line 4 lt 5 lc 3 lw 2
set style line 5 lt 5 lc 5 lw 2
set style line 6 lt 5 lc 2 lw 2
set style line 7 lt 5 lc 6 lw 2
set style line 8 lt 5 lc 8 lw 2
set style line 9 lt 5 lc 1 lw 2
set style line 10 lt 5 lc 4 lw 2
set style line 11 lt 5 lc 7 lw 2
set style line 12 lt 5 lc 9 lw 2
set style line 13 pt 6 ps 1 lc 8
set style line 14 pt 6 ps 1 lc 2
set style line 15 pt 6 ps 1 lc 4
set style line 16 pt 6 ps 1 lc 9

kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

MPatoHtrperbohr3 = 1e6*(0.52917720859e-10)**3/(1.602176487e-19*27.2117) # Converts MPa to Hartree/bohr^3

molperltobohr3 = 6.02214179e23*1000*(0.52917720859e-10)**3 # Converts mol/l to bohr^-3

#f(x) = m*x+b
#fit [0.00499:0.0055] f(x) 'figs/pressure-with-isotherms.dat' using 1:3 via m,b

plot [:0.0055] [-2e-6:5e-6] \
'figs/pressure-with-isotherms.dat' u 1:3 title 'T = 298K' with lines ls 4 , \
'figs/pressure-with-isotherms.dat' u 1:5 title 'T = 348K' with lines  ls 5, \
'figs/pressure-with-isotherms.dat' u 1:7 title 'T = 398K' with lines ls 6 , \
'figs/pressure-with-isotherms.dat' u 1:9 title 'T = 448K' with lines ls 7, \
'figs/pressure-with-isotherms.dat' u 1:11 title 'T = 498K' with lines ls 8, \
'figs/pressure-with-isotherms.dat' u 1:13 title 'T = 548K' with lines ls 9, \
'figs/pressure-with-isotherms.dat' u 1:15 title 'T = 598K' with lines ls 10, \
'figs/pressure-with-isotherms.dat' u 1:17 title 'T = 648K' with lines ls 11, \
'figs/pressure-with-isotherms.dat' u 1:19 title 'T = 698K' with lines ls 12, \
'figs/pressure-with-isotherms.dat' u 1:21 title 'T = 748K' with lines ls 4 , \
'figs/equation-of-state.dat' u 3:2 title 'eos' with lines ls 1 , \
'figs/equation-of-state.dat' u 4:2 notitle with lines ls 1 , \
'figs/experimental-equation-of-state.dat' u 3:2 title 'expt' with lines ls 2 , \
'figs/experimental-equation-of-state.dat' u 4:2 notitle with lines ls 2 , \
'figs/pressure-298K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 3 , \
'figs/pressure-398K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 14 , \
'figs/pressure-498K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 13 , \
'figs/pressure-598K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 15 , \
'figs/pressure-698K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 16 
#m*x+b title 'fit at T = 298K' with lines

set size 0.36,0.32        # The second one (inset)
set origin 0.1,0.65
unset xlabel
unset ylabel
set xtics 4e-5 font "Helvetica,12"
set ytics 1e-7 font "Helvetica,12"

plot [:8e-5] [:] \
'figs/pressure-with-isotherms.dat' u 1:3 notitle with lines ls 4 , \
'figs/pressure-with-isotherms.dat' u 1:5 notitle with lines  ls 5, \
'figs/pressure-with-isotherms.dat' u 1:7 notitle with lines ls 6 , \
'figs/pressure-with-isotherms.dat' u 1:9 notitle with lines ls 7, \
'figs/pressure-with-isotherms.dat' u 1:11 notitle with lines ls 8, \
'figs/pressure-with-isotherms.dat' u 1:13 notitle with lines ls 9, \
'figs/pressure-with-isotherms.dat' u 1:15 notitle with lines ls 10, \
'figs/pressure-with-isotherms.dat' u 1:17 notitle with lines ls 11, \
'figs/pressure-with-isotherms.dat' u 1:19 notitle with lines ls 12, \
'figs/pressure-with-isotherms.dat' u 1:21 notitle with lines ls 4 , \
'figs/equation-of-state.dat' u 3:2 notitle with lines ls 1 , \
'figs/equation-of-state.dat' u 4:2 notitle with lines ls 1 , \
'figs/experimental-equation-of-state.dat' u 3:2 notitle with lines ls 2 , \
'figs/experimental-equation-of-state.dat' u 4:2 notitle with lines ls 2 , \
'figs/pressure-298K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 3 , \
'figs/pressure-398K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 14 , \
'figs/pressure-498K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 13 , \
'figs/pressure-598K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 15 , \
'figs/pressure-698K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 16     