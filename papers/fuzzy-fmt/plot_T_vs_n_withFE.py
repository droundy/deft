#!/usr/bin/python2
#This program is for plotting Temperature vs Density
#Run this program from /deft/papers/fuzzy-fmt by entering (example):
# ./plot_T_vs_n_withFE.py [filename.dat]

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a T vs n scatter plot.')

parser.add_argument('filedat', metavar='datafile', type=str,
                    help='file with data to plot') 
args=parser.parse_args()

thisdata = np.loadtxt(args.filedat)


density=thisdata[:,1]     
T=thisdata[:,0]
FE=thisdata[:,2] 

for n in [0.9, 0.96, 1.0, 1.1]:
     plt.plot(T[density == n], FE[density == n], '-', label="n={}".format(n))
plt.legend(loc='best')
#plt.scatter(x_axis, y_axis, c=z_axis)
#levels=[-1,0, 10]
#plt.colorbar()
#plt.title(plot_title)
plt.xlabel('T')
plt.ylabel('Free energy difference')
#plt.savefig(plot_name)

plt.show()
