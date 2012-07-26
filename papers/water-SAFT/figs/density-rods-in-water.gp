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

set terminal postscript eps enhanced color "Helvetica" 20 size 4.2,4
set output 'figs/density-rods-in-water.eps'

reset
unset arrow
set view map

#set key inside bottom
#set title 'Two hydrophobic rods - Density'

set multiplot 
#set pm3d map
#set palette color positive
#set ticslevel 0
#set samples 50; set isosamples 50
#set palette rgbformulae 22,13,-31
set palette defined ( 0 "white", 0.4 "#aaffaa", 0.85 "cyan" , 1.10 "blue", 1.2 "black" )

set size 0.95,0.58          # The bottom plot
set origin 0,0
set xlabel 'y (nm)'
set ylabel 'z (nm)'
set rmargin 4
set lmargin 3
set tmargin 0
set bmargin 2
set colorbox user origin 0.81,0.128 size 0.04,0.769
set cblabel 'Density (g/mL)'
set cbrange [0:1.2]
set cbtics 0.1
set title  'After transition'

d = 0.6
set style line 1 lt 1 lw 2 lc rgb "white"
set style arrow 1 heads filled size screen 0.015, 15, 5 front ls 1
set arrow from -d/2, 0, 1 to d/2, 0, 1 as 1
set label "d" at -0.07, 0.2 rotate by 0 font 'Helvetica, 20' textcolor rgb "white" front

set style fill empty border 7
set obj 2 circle at (d/2+.5), 0   size 0.5 lw 3 front
set obj 3 circle at -(d/2+0.5), 0   size 0.5 lw 3 front

gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density
nm=18.8972613 # 1 nm in atomic units

splot [-2.4:2.4] [-1.2:1.2] [:] \
'figs/rods-slice-01.0-00.6.dat' u ($2/nm):($3/nm):($4/gpermL) notitle with pm3d

unset arrow

set size 0.95,0.53          # The top plot
set origin 0,0.48
unset xlabel
unset xtics
set ylabel 'z (nm)'
set rmargin 4
set lmargin 3
set tmargin -2
set bmargin 0
unset colorbox

d = 0.5
set style line 2 lt 1 lw 2 lc rgb "#black"
set style arrow 2 heads filled size screen 0.015, 15, 5 front ls 2
set arrow from -d/2, 0, 1 to d/2, 0, 1 as 2
set label "d" at -0.07, 0.2 rotate by 0 font 'Helvetica, 20' front

set style fill empty border 7
set obj 2 circle at (d/2+.5), 0   size 0.5 lw 3 front
set obj 3 circle at -(d/2+0.5), 0   size 0.5 lw 3 front

set title  'Before transition'

splot [-2.4:2.4] [-1.2:1.2] [:] \
'figs/rods-slice-01.0-00.5.dat' u ($2/nm):($3/nm):($4/gpermL) notitle with pm3d
