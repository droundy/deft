#!/usr/bin/python2

#This program creates plots from *best.dat or *best_tensor.dat data files
#which are generated as output data files by figs/new-melting.cpp
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-isotherm.py --kT [temp] --tensor(optional)

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser(description='Plots gw vs n and FEdiff vs n.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--tensor', action='store_true',
                    help='--tensor for use tensor weight')

args=parser.parse_args()

kT=args.kT

n = []
gw = []
fe_difference = []

if args.tensor :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best_tensor.dat' % kT)))    
else :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
for f in files:
    data = np.loadtxt(f)
    n.append(data[1])
    gw.append(data[3])
    fe_difference.append(data[6])

plt.figure()
plt.plot(n, gw, '.-')
plt.xlabel('n')
plt.ylabel('gw')

plt.figure()
plt.axhspan(0.0, -0.2, color='lightblue', alpha=0.15, lw=0)
plt.title("kT=%g" % (kT))
plt.plot(n, fe_difference, '.-')
plt.axhline(0, color='k')
plt.xlabel('n')
plt.ylabel('FEdiff')

plt.show()
