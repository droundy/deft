#!/usr/bin/python

from __future__ import division
import matplotlib, sys

if len(sys.argv) < 2 or sys.argv[1]!="show" :
  matplotlib.use('Agg')
import pylab, numpy, random
from pylab import *
#import Scientific.Functions.LeastSquares as ls

from scipy.optimize import leastsq

ff = pylab.array([.05, .1, .15, .2, .25, .3, .35, .4, .45, .5])
slope = pylab.zeros(len(ff))
hsigma = pylab.zeros(len(ff))

def read_ghs(base, ff):
    mcdatafilename = "%s-0.%02d.dat" % (base, 100*ff)
    mcdata = numpy.loadtxt(mcdatafilename)
    print 'Using', mcdatafilename, 'for filling fraction', ff
    r_mc = mcdata[:,0]
    n_mc = mcdata[:,1]
    ghs = n_mc/ff
    return r_mc, ghs

for i in range(len(ff)):
  print 'working with', i
  r,ghs = read_ghs('figs/gr', ff[i])
  hsigma[i] = (1-ff[i]/2)/(1-ff[i])**3 - 1
  skipnum = 3
  slope[i] = (ghs[skipnum] - ghs[0])/(r[skipnum] - r[0])
  plot(r, ghs[0] + (r-2)*slope[i])
  plot(r, ghs)
pylab.xlim(2,2.2)
pylab.ylim(0, 6)

pylab.figure()
pylab.plot(hsigma, slope, 'b-')
pylab.plot(hsigma, -hsigma-hsigma**2, 'r-')
axhline(y=0)
pylab.legend(loc='best')
pylab.show()
