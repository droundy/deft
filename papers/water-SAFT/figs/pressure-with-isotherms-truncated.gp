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
set output 'figs/pressure-with-isotherms-truncated.eps'

unset key
# set key noauto outside bottom
# set key font "Helvetica,14"
# set key horizontal

set multiplot

set size 1,1          # The first plot (host plot)
set origin 0,0
set xlabel 'Density (g/mL)'
set ylabel 'Pressure (atm)'

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

set style line 13 pt 6 ps 1 lw 0 lt 0 lc 8
set style line 14 pt 6 ps 1 lw 0 lt 0 lc 2
set style line 15 pt 6 ps 1 lw 0 lt 0 lc 4
set style line 16 pt 6 ps 1 lw 0 lt 0 lc 9
set style line 17 pt 6 ps 1 lw 0 lt 0 lc 6
set style line 18 pt 6 ps 1 lw 0 lt 0 lc 5

kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin

MPatoHtrperbohr3 = 1e6*(0.52917720859e-10)**3/(1.602176487e-19*27.2117) # Converts MPa to Hartree/bohr^3

molperltobohr3 = 6.02214179e23*1000*(0.52917720859e-10)**3 # Converts mol/l to bohr^-3

atm=3.4439674e-09 # in atomic units
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

#f(x) = m*x+b
#fit [0.00499:0.0055] f(x) 'figs/pressure-with-isotherms.dat' using 1:3 via m,b

set label "298K" at 0.985,10 rotate by 89 font 'Helvetica,14'
set label "348K" at 0.945,10 rotate by 89 font 'Helvetica,14'
set label "398K" at 0.905,15 rotate by 89 font 'Helvetica,14'
set label "448K" at 0.860,30 rotate by 89 font 'Helvetica,14'
set label "498K" at 0.805,50 rotate by 89 font 'Helvetica,14'
set label "548K" at 0.745,10 rotate by 88 font 'Helvetica,14'
set label "598K" at 0.665, 30 rotate by 86 font 'Helvetica,14'
set label "648K" at 0.535, 110 rotate by 70 font 'Helvetica,14'
set label "698K" at 0.24, 360 rotate by 10 font 'Helvetica,14'
set label "748K" at 0.175, 400 rotate by 55 font 'Helvetica,14'
#set arrow from 0.001,1e-6 to 0.002,1e-6

set logscale y

plot [:1.05] [0.01:150] \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($3/atm) title 'T = 298K' with lines ls 4 , \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($7/atm) title 'T = 398K' with lines ls 6 , \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($11/atm) title 'T = 498K' with lines ls 8, \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($15/atm) title 'T = 598K' with lines ls 10, \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($5/atm) title 'T = 348K' with lines  ls 5, \
'figs/pressure-with-isotherms.dat' u ($1/gpermL):($9/atm) title 'T = 448K' with lines ls 7, \
'figs/equation-of-state.dat' u ($3/gpermL):($2/atm) title 'eos' with lines ls 1 , \
'figs/equation-of-state.dat' u ($4/gpermL):($2/atm) notitle with lines ls 1 , \
'figs/experimental-equation-of-state.dat' u ($3/gpermL):($2/atm) title 'expt' with lines ls 2 , \
'figs/experimental-equation-of-state.dat' u ($4/gpermL):($2/atm) notitle with lines ls 2 , \
'figs/pressure-298K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 3 , \
'figs/pressure-348K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 18 , \
'figs/pressure-398K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 14 , \
'figs/pressure-448K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 17 , \
'figs/pressure-498K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 13 , \
'figs/pressure-598K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 15# , \
#'figs/pressure-698K.dat' u ($3*molperltobohr3/gpermL):($2*MPatoHtrperbohr3/atm) notitle with lp ls 16 
#'figs/pressure-with-isotherms.dat' u ($1/gpermL):($19/atm) title 'T = 698K' with lines ls 12, \
#'figs/pressure-with-isotherms.dat' u ($1/gpermL):($21/atm) title 'T = 748K' with lines ls 4 , \


#'figs/pressure-with-isotherms.dat' u ($1/gpermL):($13/atm) title 'T = 548K' with lines ls 9, \
#'figs/pressure-with-isotherms.dat' u ($1/gpermL):($17/atm) title 'T = 648K' with lines ls 11, \

#m*x+b title 'fit at T = 298K' with lines

# set size 0.36,0.32        # The second one (inset)
# set origin 0.1,0.65
# unset xlabel
# unset ylabel
# set xtics 4e-5 font "Helvetica,12"
# set ytics 1e-7 font "Helvetica,12"

# plot [:8e-5] [:] \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):3 notitle with lines ls 4 , \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):5 notitle with lines  ls 5, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):7 notitle with lines ls 6 , \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):9 notitle with lines ls 7, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):11 notitle with lines ls 8, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):13 notitle with lines ls 9, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):15 notitle with lines ls 10, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):17 notitle with lines ls 11, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):19 notitle with lines ls 12, \
# 'figs/pressure-with-isotherms.dat' u ($1/gpermL):21 notitle with lines ls 4 , \
# 'figs/equation-of-state.dat' u 3:2 notitle with lines ls 1 , \
# 'figs/equation-of-state.dat' u 4:2 notitle with lines ls 1 , \
# 'figs/experimental-equation-of-state.dat' u 3:2 notitle with lines ls 2 , \
# 'figs/experimental-equation-of-state.dat' u 4:2 notitle with lines ls 2 , \
# 'figs/pressure-298K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 3 , \
# 'figs/pressure-398K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 14 , \
# 'figs/pressure-498K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 13 , \
# 'figs/pressure-598K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 15 , \
# 'figs/pressure-698K.dat' u ($3*molperltobohr3):($2*MPatoHtrperbohr3) notitle with points ls 16
