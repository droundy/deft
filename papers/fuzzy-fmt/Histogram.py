#!/usr/bin/python2
#Run from /deft/papers/fuzzy-fmt by entering ./Histogram.py

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

data_file="histogram_relerror.dat"
thisdata = np.loadtxt(data_file)

#x=thisdata[:,6] 

x=[1.50E-03,
1.23E-03,
1.12E-03,
8.06E-04,
6.12E-04,
6.00E-04,
5.57E-04,
5.61E-04,
5.40E-04,
5.08E-04,
4.66E-04,
4.12E-04,
4.00E-04,
3.86E-04,
3.24E-04,
2.20E-04,
2.19E-04,
1.51E-04,
1.35E-04,
7.44E-05,
6.77E-05,
1.54E-05,
-4.22E-05,
-5.87E-05,
-6.15E-05,
-1.59E-04,
-1.61E-04,
-2.39E-04,
-2.40E-04,
-2.57E-04,
-2.60E-04,
-2.81E-04,
-3.02E-04,
3.24E-04,
-3.70E-04,
-3.70E-04,
-3.84E-04,
-4.60E-04,
-4.80E-04,
-5.21E-04,
-5.49E-04,
-5.93E-04,
-6.15E-04,
-6.39E-04,
-7.44E-04,
-9.03E-04,
-1.06E-03,
-1.66E-03,
-1.70E-03]

num_bins=6
#num_bins=100

n,bins,patches=plt.hist(x,num_bins,facecolor='blue', alpha=0.5)
x_label="relerror"
plt.xlabel(x_label)
plt.show()
