#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./nm_plot_histdata.py [filename.dat]


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")
parser.add_argument('filedat', metavar='datafile', type=str, nargs='*',
                    help='file with data to plot') 
args=parser.parse_args()

markers = {
   0.1: '^',
   0.5: 'o',
   0.01: '+',
   0.05: 'v',
}

mc_constant = {}
percent_neg = {}
estimated_uncertainty = {}
mean_error = {}
error_width = {}
mean_error_uncertainty = {}

for f in args.filedat:
   print 'Reading file', f
   thisdata = np.loadtxt(f)

   gw = thisdata[0,12]
   mc_prefactor = thisdata[0,11]
   if gw not in mc_constant.keys():
      mc_constant[gw] = []
      percent_neg[gw] = []
      estimated_uncertainty[gw] = []
      mean_error[gw] = []
      error_width[gw] = []
      mean_error_uncertainty[gw] = []

   mc_constant[gw].append(thisdata[0,10])
   rel_error=thisdata[:,9]
   num_seeds = len(thisdata)

   percent_neg[gw].append(100.0*len(rel_error[rel_error<0]) / float(len(rel_error)))
   estimated_uncertainty[gw].append(100.0*len(rel_error[rel_error<0]) / float(len(rel_error))/np.sqrt(len(rel_error)))

   mean_error[gw].append(rel_error.mean())
   error_width[gw].append(rel_error.std())
   mean_error_uncertainty[gw].append(rel_error.std()/np.sqrt(len(rel_error)))

for gw in mc_constant.keys():
   plt.errorbar(mc_constant[gw], percent_neg[gw], estimated_uncertainty[gw], fmt=markers[gw], label='gw=%g' % gw)
plt.legend()

plt.axhspan(60, 40, color='green', alpha=0.15, lw=0)
plt.axhspan(49.5, 50.5, color='black', alpha=0.15, lw=0)

plt.title('Histogram results with mc_prefactor %g' % mc_prefactor)
plt.ylabel('percent of negative relerrors')
plt.xlabel('mc_constant')

plt.figure()

plt.axhspan(-1e-4, 1e-4, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='green')

for gw in mc_constant.keys():
   plt.errorbar(mc_constant[gw], mean_error[gw], mean_error_uncertainty[gw], fmt=markers[gw], label='gw=%g' % gw)
plt.legend()

plt.title('Histogram results with mc_prefactor %g' % mc_prefactor)
plt.ylabel('mean relative error')
plt.xlabel('mc_constant')

plt.show()
