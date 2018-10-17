#!/usr/bin/python2
#Run from /deft/papers/fuzzy-fmt by entering ./Histogram.py [filename.dat]

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

parser.add_argument('filedat', metavar='histogram datafile', type=str,
                    help='data file used to create histogram') 

args=parser.parse_args()

thisdata = np.loadtxt(args.filedat)

mc_constant = thisdata[0,10]
gw = thisdata[0,12]
rel_error=thisdata[:,9]
num_seeds = len(thisdata)

error_range = max(rel_error.max(), -rel_error.min())
_,bins,_ = plt.hist(rel_error,range=(-error_range, error_range), facecolor='blue')
plt.hist(rel_error[rel_error<0],bins,facecolor='green')
plt.title(r'%.1f%% $\pm$ %.1f%% negative  %s' %
          (100.0*len(rel_error[rel_error<0]) / float(len(rel_error)),
           100.0*len(rel_error[rel_error<0]) / float(len(rel_error))/np.sqrt(len(rel_error)), args.filedat))
percent_neg=100.0*len(rel_error[rel_error<0]) / float(len(rel_error))
print('%s %g %g %g    #%.1f%% +- %.1f%% negative  %s' %
          (percent_neg, gw, mc_constant, num_seeds, 100.0*len(rel_error[rel_error<0]) / float(len(rel_error)),
           100.0*len(rel_error[rel_error<0]) / float(len(rel_error))/np.sqrt(len(rel_error)), args.filedat))
#print(bins)
plt.xlabel("relerror")
plt.savefig(args.filedat+".png" )
plt.show()

