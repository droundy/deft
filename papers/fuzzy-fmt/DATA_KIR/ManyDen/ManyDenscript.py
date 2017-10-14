#!/usr/bin/python2

import os

os.system('rm newmeltdataout.dat')
os.system('rm plot.dat')

#firstlongrun for rdensity in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.9]:
#run2 for rdensity in [1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78]:
for rdensity in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]:
#for rdensity in [0.9, 1.0, 1.2] :
    os.system('figs/new-melting.mkdat 2 %g -1 -1' % (rdensity)) 
    os.system('cat newmeltbestdata.dat >> plot.dat')    

os.system('gnuplot newmeltbestdata_plot.gnu')
