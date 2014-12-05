#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy.special import erf
import matplotlib.lines as mlines
import os

import styles

if len(sys.argv) != 3:
    print("Usage:  " + sys.argv[0] + ' filling-fraction temperatures')
    exit(1)

reduced_density = float(sys.argv[1])


all_temperatures = eval(sys.argv[2])


print 'all_temperatures are', all_temperatures

#Constants and variables
dT = .001
V0 = 1
R = 1# in Angstroms
density = arange(0, .8 - .001/2, .001)/(4*pi/3)
eta = density*4*pi/3
P_cs = density*(1+eta+eta**2)/(1-eta)**3


print(reduced_density)
density = (reduced_density/100)/(4*pi/3)
phs = density*(1+(reduced_density/100)+(reduced_density/100)**2)/(1-(reduced_density/100))**3

mysymbol_names = []
mysymbol_lines = []

mylines = []


#for temp in all_temperatures:
  # input: ['figs/mc_testp_wca-0.%02d00-%.4f.dat.gradial' % (reduced_density, temp) for temp in all_temperatures]
temp = all_temperatures
fname = 'figs/mc_testp_wca-%.4f-%.4f.dat.gradial' % (reduced_density,temp)
if os.path.exists(fname):
  print 'found', fname
  g = loadtxt(fname)
  line = plot(g[:,0], g[:,1])
  xlim(xmax=8)
else:
  print 'could not find', fname

title('Radial distribution function at $n^* = %g$' % (reduced_density/100))
xlabel('radius')
ylabel('g')
blue_line = mlines.Line2D([], [], color='blue', marker='*',
                          markersize=15, label='Blue stars')

savefig('figs/testp-radial-distribution-%02d.pdf' % (reduced_density), bbox_inches=0)
print('figs/testp-radial-distribution-%02d.pdf' % (reduced_density))

show()
