#!/usr/bin/python2
#Run from /deft/papers/fuzzy-fmt by entering ./Histogram.py [filename.dat] --gw [] --mcconstant [] --seeds []

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

parser.add_argument('filedat', metavar='histogram datafile', type=str,
                    help='data file used to create histogram') 
parser.add_argument('--gw', metavar='gwidth', type=float,
                    help='gwidth') 
parser.add_argument('--mcconstant', metavar='mcconstant', type=float,
                    help='mc_constant') 
parser.add_argument('--seeds', metavar='seeds', type=float,
                    help='number of seeds ran') 
                                        
args=parser.parse_args()

thisdata = np.loadtxt(args.filedat)

x=thisdata[:,9]

error_range = max(x.max(), -x.min())
_,bins,_ = plt.hist(x,range=(-error_range, error_range), facecolor='blue')
plt.hist(x[x<0],bins,facecolor='green')
plt.title(r'%.1f%% $\pm$ %.1f%% negative  %s' %
          (100.0*len(x[x<0]) / float(len(x)),
           100.0*len(x[x<0]) / float(len(x))/np.sqrt(len(x)), args.filedat))
percent_neg=100.0*len(x[x<0]) / float(len(x))
print('%s %g %g %g    #%.1f%% +- %.1f%% negative  %s' %
          (percent_neg, args.gw, args.mcconstant, args.seeds, 100.0*len(x[x<0]) / float(len(x)),
           100.0*len(x[x<0]) / float(len(x))/np.sqrt(len(x)), args.filedat))
x_label="relerror"
#print(bins)
plt.xlabel(x_label)
plt.savefig(args.filedat+".png" )
plt.show()

