#!/usr/bin/python2

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.interpolate

parser = argparse.ArgumentParser()
parser.add_argument('directory', metavar='exisiting_data_directory', type=str,
                    help='exisiting directory for data files')
parser.add_argument('--temp', help='the temperature')
parser.add_argument('--density', help='the density')
parser.add_argument('--de-max', help='maximum energy difference to plot')
parser.add_argument('--de-min', help='minimum energy difference to plot')

args=parser.parse_args()

data_directory=args.directory

data_file=data_directory+"/plot_alldat.dat"
thisdata = np.loadtxt(data_file)
print
print "Using data from file:"+data_file
print

T = thisdata[:,0]
fv = thisdata[:,2]
gw = thisdata[:,3]

free_energy_difference = thisdata[:,6]

#plt.show()
#fig, ax = plt.subplots()
#cmap = mp1.colors.ListedColormap(['red', 'blue'])
#bounds = [-2, 0, 2]
#norm =mp1.colors.BoundaryNorm(bounds, cmap.N)
#cb2 = mp1.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, boundaries=[0] + bounds + [13], extend='both', ticks=bounds, spacing='proportional', orientation='vertical')
#cb2.set_label('DiffFE')
#fig.show
##contour=plt.contour(x,y,zi, levels, colors='k')
##plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize= 12)
#contour_filled = plt.contourf(x,y,zi, levels)
#plt.colorbar()

#plt.show()
plt.title("Difference in Free Energy, kT=2, n=0.92")
plt.xlabel("fv")
plt.ylabel("gw")

plt.scatter(fv[okay_data],gw[okay_data],c=free_energy_difference[okay_data])
plt.colorbar()
plt.show()
