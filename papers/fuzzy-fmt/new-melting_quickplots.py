#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./new-melting_quickplots.py  --kT value --n value  --d directory
#Generates 8 plot.dat files and plots

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

args=parser.parse_args()
directory=args.d
kT=args.kT
n=args.n

#Try values of n: 0.95, 1.08, 1.2  and  kT: 1, 2, 20, 200

for gwidth in [0.001, 0.01, 0.1, .2]: 
    os.system('figs/new-melting.mkdat --kT %g --n %g --d %s_kT%g_n%g_gw%g --fvstart 0 --fvend .9 --fvstep 0.01 --gw %g  --dx 0.01' % (kT, n, directory, kT, n, gwidth, gwidth))
    os.system('cat %s_kT%g_n%g_gw%g/*alldat.dat >> %s_kT%g_n%g_gw%g/plot.dat' % (directory, kT, n, gwidth, directory, kT, n, gwidth))
    os.system('./new-melting_anyplot_script.py %s_kT%g_n%g_gw%g --ftemp %g --ydiff --xfv --ptname n%g_gw%g' % (directory, kT, n, gwidth, kT, n, gwidth))
    os.system('cp %s_kT%g_n%g_gw%g/plot_DiffFEvsfv_kT%g_n%g_gw%g.png %s_plots/plot_DiffFEvsfv_kT%g_n%g_gw%g.png' % (directory, kT, n, gwidth, kT, n, gwidth, directory, kT, n, gwidth))

for fv in [0, 0.01, 0.1, 0.2]:
    os.system('figs/new-melting.mkdat --kT %g --n %g --d %s_kT%g_n%g_fv%g --fv %g --gwstart .001 --gwend .3 --gwstep .005  --dx 0.01' % (kT, n, directory, kT, n, fv, fv))
    os.system('cat %s_kT%g_n%g_fv%g/*alldat.dat >> %s_kT%g_n%g_fv%g/plot.dat' % (directory, kT, n, fv, directory, kT, n, fv))
    os.system('./new-melting_anyplot_script.py %s_kT%g_n%g_fv%g --ftemp %g --ydiff --xgw --ptname n%g_fv%g' % (directory, kT, n, fv, kT, n, fv))
    os.system('cp %s_kT%g_n%g_fv%g/plot_DiffFEvsgw_kT%g_n%g_fv%g.png %s_plots/plot_DiffFEvsgw_kT%g_n%g_fv%g.png' % (directory, kT, n, fv, kT, n, fv, directory, kT, n, fv))

