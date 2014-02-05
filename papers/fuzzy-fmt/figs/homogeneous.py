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

#Constants and variables
#k_b = 8.6173324*10**(-5) # in eV
dT = .001
#Temp_max = 600 #in Kelvin
#Temp = arange(.001, .1 + dT/2, dT)
V0 = 1
R = 1# in Angstroms
density = arange(0, .8 - .001/2, .0001)/(4*pi/3)

eta = density*4*pi/3
#P_cs = density*.001*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
#plot(eta,P_cs/.001, 'k',linewidth=2, label = 'Hard spheres')


colors = { 0.1: 'r', 0.01: 'm', 0.001: 'b', 0.0001: 'c', 0.00001: 'g' }

Temp = 0.00001
eta = density*4*pi/3
#P_cs = density*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
P_cs = density
while Temp <= .011:
  betaV0 = V0/Temp

  gamma = betaV0*(4+sqrt((4+4*sqrt(pi/betaV0)-2*pi*sqrt(pi/betaV0)+pi/betaV0)**2-4*pi**3/betaV0)+4*sqrt(pi/betaV0)-2*pi*sqrt(pi/betaV0)+pi/betaV0)/(2*pi**2)
  sg = sqrt(gamma)

  #Integrals for the different weighted densities
  W3 = (-pi*R/(3*gamma**(3/2)*(sqrt(pi*gamma) -1)))*(2*sg*(8*(1 + gamma) - exp(-gamma)*(2*gamma+5))-sqrt(pi)*(4*gamma**2+12*gamma+3)*erf(sg))
  W2 = (2*pi*R**2/(gamma*(sqrt(pi*gamma)-1)))*(sqrt(pi*gamma)*(2*gamma + 3)*erf(sg) - 6*gamma+2*exp(-gamma)*(gamma + 1) -  2)
  W1 = R/(2*sg*(sqrt(pi*gamma)-1))*(sqrt(pi)*(2*gamma + 1)*erf(sg) + sg*(2*exp(-gamma)-4))
  W0 = (sg/(sqrt(gamma*pi)-1))*(sqrt(pi)*erf(sg) + (exp(-gamma)-1)/sg)

  #Weighted densities
  n0 = density*W0
  n1 = density*W1
  n2 = density*W2
  n3 = density*W3

  Phi_1 = -n0*log(1-n3)
  Phi_2 = n1*n2/(1-n3)
  Phi_3 = n2**3/3/(8*pi*(1-n3)**2)

  Phi = Phi_1 + Phi_2 + Phi_3

  alpha = 1-n3
  dPhi1_dn = Phi_1/density + density*W0*W3/alpha
  dPhi2_dn = (W1*n2+W2*n1)/alpha + n1*n2/alpha**2*W3
  dPhi3_dn = n2**2/(8*pi*alpha**2)*W2 + n2**3/12/pi/alpha**3*W3
  dPhi_dn = dPhi1_dn + dPhi2_dn + dPhi3_dn


  pressure = Temp*(density*dPhi_dn - Phi) + density*Temp 
  plot(density*(4*pi/3),pressure/P_cs/Temp, colors[Temp]+'-', label = 'T/V0=%g' %Temp)
  Temp = Temp*10 # sqrt(10)

erfdata = loadtxt('figs/homogeneous.dat')
erftemps = arange(0.01, 0.03, 0.01)
erfeta = erfdata[:,0]
erf_density = erfeta/(4*pi/3)
for j in arange(1,len(erfdata[0,:])):
  erfpressure = erfdata[:,j]
  plot(erfeta,erfpressure, '--', label='kT=%g' % (0.01*(j-1)))
  

for ff in arange(0.1,0.81, 0.1):
  density = ff/(4*pi/3)
  phs = density*(1+ff+ff**2)/(1-ff)**3
  phs = density
  for temp in [0.1, 0.01, 0.001, 0.0001]:
    fname = 'figs/mc-%.4f-%.4f.dat.prs' % (ff, temp)
    if os.path.exists(fname):
      print 'found', fname
      p = loadtxt(fname)
      plot(ff, p/temp/phs, colors[temp] + 'o')
    else:
      print 'could not find', fname

#plot(density*(4*pi/3), density, label = 'ideal gas')
erfdata = loadtxt('figs/homogeneous.dat')
ylim(1,15.0)
xlim(xmax=0.51)
#mcdata = loadtxt('figs/mc-soft-homogenous-20-382-1.00000.dat.prs')
#plot(mcdata[:,1],mcdata[:,0],'*')
xlabel('packing fraction')
ylabel('pressure / hard-sphere pressure')
legend(loc = 'best')
savefig('figs/p-vs-packing.pdf', bbox_inches=0)
show()
