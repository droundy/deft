#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./new-melting_quickplots.py  --kT 2 --n 1.2  --d carrots --dx 0.001 (optional)
#after making directory (ie. carrots_plots)
#Generates 8 plot.dat files and plots and puts them in the existing directory carrots_plots

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

parser.add_argument('--d', metavar='exisiting_data_directory', type=str,
                    help='exisiting directory for data files') 
parser.add_argument('--kT', metavar='temperature', type=float,
                    help='create plot at this temperature')
parser.add_argument('--n', metavar='reduced_density', type=float,
                    help='create plot at this reduced density')
parser.add_argument('--dx', metavar='dx', type=float,
                    help='dx')

args=parser.parse_args()
directory=args.d
kT=args.kT
n=args.n

dx=0.01
if args.dx:
    dx=args.dx

#Try values of n: 0.95, 1.08, 1.2  and  kT: 1, 2, 20, 200

for gwidth in [0.001, 0.01, 0.1, .2]: 
    os.system('figs/new-melting.mkdat --kT %g --n %g --d %s_kT%g_n%g_gw%g_dx%g --fvstart 0 --fvend .1 --fvstep 0.01 --gw %g  --dx %g' % (kT, n, directory, kT, n, gwidth, dx, gwidth, dx))
    os.system('cat %s_kT%g_n%g_gw%g_dx%g/*alldat.dat >> %s_kT%g_n%g_gw%g_dx%g/plot.dat' % (directory, kT, n, gwidth, dx, directory, kT, n, gwidth, dx))
    os.system('./new-melting_anyplot_script.py %s_kT%g_n%g_gw%g_dx%g --ftemp %g --ydiff --xfv --ptname n%g_gw%g_dx%g' % (directory, kT, n, gwidth, dx, kT, n, gwidth, dx))
    os.system('cp %s_kT%g_n%g_gw%g_dx%g/plot_DiffFEvsfv_kT%g_n%g_gw%g_dx%g.png %s_plots/plot_DiffFEvsfv_kT%g_n%g_gw%g_dx%g.png' % (directory, kT, n, gwidth, dx, kT, n, gwidth, dx, directory, kT, n, gwidth, dx))

for fv in [0, 0.01, 0.1]:
    os.system('figs/new-melting.mkdat --kT %g --n %g --d %s_kT%g_n%g_fv%g_dx%g --fv %g --gwstart .001 --gwend .7 --gwstep .01  --dx %g' % (kT, n, directory, kT, n, fv, dx, fv, dx))
    os.system('cat %s_kT%g_n%g_fv%g_dx%g/*alldat.dat >> %s_kT%g_n%g_fv%g_dx%g/plot.dat' % (directory, kT, n, fv, dx, directory, kT, n, fv, dx))
    os.system('cat %s_kT%g_n%g_fv%g_dx%g/*alldat.dat >> %s_kT%g_n%g_fv%g_dx%g/plot.dat' % (directory, kT, n, fv, dx, directory, kT, n, fv, dx))
    os.system('./new-melting_anyplot_script.py %s_kT%g_n%g_fv%g_dx%g --ftemp %g --ydiff --xgw --ptname n%g_fv%g_dx%g' % (directory, kT, n, fv, dx, kT, n, fv, dx))
    os.system('cp %s_kT%g_n%g_fv%g_dx%g/plot_DiffFEvsgw_kT%g_n%g_fv%g_dx%g.png %s_plots/plot_DiffFEvsgw_kT%g_n%g_fv%g_dx%g.png' % (directory, kT, n, fv, dx, kT, n, fv, dx, directory, kT, n, fv, dx))

