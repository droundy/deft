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
#arg reduced_density = [10, 20, 30, 40, 50]

all_temperatures = eval(sys.argv[2])
#arg all_temperatures = [[1.0, 0.1, 0.01, 0.001]]

print 'all_temperatures are', all_temperatures

#Constants and variables
#k_b = 8.6173324*10**(-5) # in eV
dT = .001
#Temp_max = 600 #in Kelvin
#Temp = arange(.001, .1 + dT/2, dT)
V0 = 1
#betaV0 = V0/Temp
R = 1# in Angstroms
density = arange(0, .8 - .001/2, .001)/(4*pi/3)
#gamma = 2*((sqrt(pi*betaV0)+sqrt(pi*betaV0-16*sqrt(betaV0)))/8)**2
#sg = sqrt(gamma)

eta = density*4*pi/3
#P_cs = density*.001*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
#plot(eta,P_cs/.001, 'k',linewidth=2, label = 'Hard spheres')


density = (reduced_density/100)/(4*pi/3)
phs = density*(1+(reduced_density/100)+(reduced_density/100)**2)/(1-(reduced_density/100))**3

mysymbol_names = []
mysymbol_lines = []

mylines = []

for temp in all_temperatures:
  # input: ['figs/radial-wca-%06.4f-%04.2f.dat' % (temp, reduced_density/100.0) for temp in all_temperatures]
  fname = 'figs/radial-wca-%06.4f-%04.2f.dat' % (temp, reduced_density/100.0)
  data = loadtxt(fname)
  r = data[:,0]
  nreduced_density = data[:,1]
  g = nreduced_density/(reduced_density/100.0)
  line = plot(r, g, styles.dftwca[temp])#, label='WCA DFT $kT/V_{max}$ = %g' % temp)
  mylines += line
  if temp == all_temperatures[-1]:
      mysymbol_lines += line
      mysymbol_names += ['WCA DFT']
  #xlim(xmax=floor(max(g[:,0])))

for temp in all_temperatures:
  # input: ['figs/mcwca-0.%02d00-%.4f.dat.gradial' % (reduced_density, temp) for temp in all_temperatures]
  fname = 'figs/mcwca-0.%02d00-%.4f.dat.gradial' % (reduced_density, temp)
  if os.path.exists(fname):
    print 'found', fname
    g = loadtxt(fname)
    line = plot(g[:,0], g[:,1], styles.mcwca[temp])#, label = 'WCA MC $kT/?$ = %g' % temp)
    #xlim(xmax=floor(max(g[:,0])))
    xlim(xmax=8)
    if temp == all_temperatures[-1]:
        mysymbol_lines += line
        mysymbol_names += ['WCA MC']
  else:
    print 'could not find', fname

title('Radial distribution function at $n^* = %g$' % (reduced_density/100))
xlabel('radius')
ylabel('g')
blue_line = mlines.Line2D([], [], color='blue', marker='*',
                          markersize=15, label='Blue stars')
legend(mylines + mysymbol_lines,
       ['$T^* = %g$' % t for t in all_temperatures] + mysymbol_names,
       loc = 'best')
savefig('figs/radial-distribution-%02d.pdf' % (reduced_density), bbox_inches=0)
print('figs/radial-distribution-%02d.pdf' % (reduced_density))

show()
