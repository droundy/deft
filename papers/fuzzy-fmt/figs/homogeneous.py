#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *
from scipy.special import erf
import os
import styles

#Constants and variables
#k_b = 8.6173324*10**(-5) # in eV
dT = .001
#Temp_max = 600 #in Kelvin
#Temp = arange(.001, .1 + dT/2, dT)
epsilon = 1
R = 1# in Angstroms
density = arange(0, .8 - .001/2, .0001)/(4*pi/3)

eta = density*4*pi/3
#P_cs = density*.001*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
#plot(eta,P_cs/.001, 'k',linewidth=2, label = 'Hard spheres')

Temp = 0.00001
eta = density*4*pi/3
#P_cs = density*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
P_cs = density

erfdata = loadtxt('figs/homogeneous.dat')
erftemp = erfdata[0,:]
erfeta = erfdata[:,0]
for j in arange(1,len(erfdata[0,:])):
  erfpressure = erfdata[:,j]
  plot(erfeta[1:]*3*2**(5.0/2.0)/4/pi, erfpressure[1:], styles.color[erftemp[j]]+'-', label='$kT/\epsilon$=%g' % erftemp[j])
#title('FIXME: Check on weighting functions for homogeneous at low $T$')

for rd in arange(0.1,1.1, 0.1):
  density = rd*2**(-5.0/2.0)
  for temp in [10.0,1.0, 0.1, 0.01, 0.001]:
    #input: 'figs/mcwca-%.4f-%.4f.dat.prs' % (rd,temp)
    fname = 'figs/mcwca-%.4f-%.4f.dat.prs' % (rd, temp)
    if os.path.exists(fname):
      print 'found', fname
      p = loadtxt(fname)
      plot(rd, p/(temp*density), styles.color[temp] + 'o')
    else:
      print 'could not find', fname

#plot(density*(4*pi/3), density, label = 'ideal gas')
ylim(ymin=1, ymax=8)
xlim(xmax=0.95)
#mcdata = loadtxt('figs/mc-soft-homogenous-20-382-1.00000.dat.prs')
#plot(mcdata[:,1],mcdata[:,0],'*')
xlabel('reduced density')
ylabel('pressure / ideal gas pressure')
legend(loc = 'best')
savefig('figs/p-vs-packing.pdf', bbox_inches=0)
show()
