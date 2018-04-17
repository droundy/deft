#!/usr/bin/python2
#This program is for plotting FE vs Density
#Run this program from /deft/papers/fuzzy-fmt by entering (example):
#  ./new-melting_FEvsDen_quickplot.py  --kT 2 --fv 0 --gw 0.001 --nstart 0.1 --nstop 2.0 --nstep 0.1 --d carrots --dx 0.001 (optional)


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

parser.add_argument('--d', metavar='exisiting_data_directory', type=str,
                    help='exisiting directory for data files') 
parser.add_argument('--kT', metavar='temperature', type=str,
                    help='create plot at this temperature')
parser.add_argument('--fv', metavar='fraction of vacancies', type=str,
                    help='create plot at this fraction of vacancies')
parser.add_argument('--gw', metavar='width of Gaussian', type=str,
                    help='create plot at this gwidth')
parser.add_argument('--nstart', metavar='start densities at', type=float,
                    help='plot densities starting at')
parser.add_argument('--nstop', metavar='stop densites at', type=float,
                    help='plot densities stopping at')
parser.add_argument('--nstep', metavar='step density by', type=float,
                    help='plot densites with step')
parser.add_argument('--dx', metavar='dx', type=str,
                    help='dx')

args=parser.parse_args()
kT=args.kT
fv=args.fv
gwidth=args.gw
nstart=args.nstart
nstop=args.nstop + args.nstep
nstep=args.nstep

dx="0.01"
if args.dx:
    dx=args.dx
    
directory=args.d+"kT"+kT+"fv"+fv+"gw"+gwidth+"dx"+dx

#for n in np.arange(nstart, nstop, nstep): 
#     os.system('figs/new-melting.mkdat --kT %s --n %s --d %s --fv %s --gw %s  --dx %s' % (kT, n, directory, fv, gwidth, dx))
     
#os.system('cat %s/*alldat.dat >> %s/plot.dat' % (directory, directory))

data_file=directory+"/plot.dat"
thisdata = np.loadtxt(data_file)

x_axis=thisdata[:,1]     
x_label="Reduced Density (n)"
x_plot="n"

y_axis=thisdata[:,6]  
y_label="Diff=(cFE-hFE)/atom" 
y_plot="DiffFE"

plot_name=directory+"/plot_"+y_plot+"vs"+x_plot+"_"+"kT"+kT+"fv"+fv+"gw"+gwidth+".png"
plot_title=y_label+" vs "+x_label+" at kT"+kT+"fv"+fv+"gw"+gwidth
a=0 
b=-50

#Plot x-axis vs y-axis
plt.axhspan(a, b, color='b', alpha=0.15, lw=0)
#plt.plot(x_axis, y_axis, color="purple")
plt.scatter(x_axis, y_axis, color="blue")
plt.title(plot_title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.savefig(plot_name)

plt.show()
