#!/usr/bin/python2
#Run from /deft/papers/fuzzy-fmt by entering ./Histogram.py

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

thisdata = np.loadtxt("histogram_relerror.dat")

x=thisdata[:,6]

error_range = max(x.max(), -x.min())
_,bins,_ = plt.hist(x,range=(-error_range, error_range), facecolor='blue')
plt.hist(x[x<0],bins,facecolor='green')
plt.title(r'%.1f%% $\pm$ %.1f%% negative' %
          (100.0*len(x[x<0]) / float(len(x)),
           100.0*len(x[x<0]) / float(len(x))/np.sqrt(len(x))))
x_label="relerror"
print(bins)
plt.xlabel(x_label)
plt.show()
