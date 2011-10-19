#!/usr/bin/gnuplot

set terminal postscript eps enhanced solid color "Helvetica" 24

set output 'vector_surface_E0_T300_x1_N16.eps'

reset
unset arrow
set xrange [-1:16]
set yrange [-1:16]
set view map
set title "Polarization Field (x=1, N=16, T=300, E=0)"
set xlabel "z"
set ylabel "y"
set palette rgbformulae 22,13,-31
splot "vector_plot_E0_T300_x1_N16.ising" u 1:2:3 with pm3d notitle






