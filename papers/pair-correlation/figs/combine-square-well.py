#!/usr/bin/python

import glob, pylab

eta = 0.3

data = pylab.loadtxt("figs/mc/a1/wallsMC-a1-pair-%02.1f-%05.3f.dat" %(eta,2.005))
z = data[:,0]
dadz = pylab.zeros_like(z)
dr = 0.01
rsquare_well = 1.79*2
for r in pylab.arange(2.005, rsquare_well, dr):
    data = pylab.loadtxt("figs/mc/a1/wallsMC-a1-pair-%02.1f-%05.3f.dat" %(eta,r))
    rmax = r + dr/2
    rmin = r - dr/2
    dadz += data[:,1]*4*pylab.pi/3*(rmax**3 - rmin**2)

data[:,1] = dadz/(4*pylab.pi*rsquare_well**3/3 - 4*pylab.pi*2.0**3/3)
pylab.savetxt("figs/mc/a1/square-well-%.1f-%.3f.dat" % (eta, rsquare_well/2), data, fmt='%.3g')
