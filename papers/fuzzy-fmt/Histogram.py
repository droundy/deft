#!/usr/bin/python2
#Run from /deft/papers/fuzzy-fmt by entering ./Histogram.py

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

data_file="histogram_relerror.dat"
thisdata = np.loadtxt(data_file)

x=thisdata[:,6]

num_bins=6
#num_bins=100

n,bins,patches=plt.hist(x,num_bins,facecolor='blue', alpha=0.5)
x_label="relerror"
plt.xlabel(x_label)
plt.show()
