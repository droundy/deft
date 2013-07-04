#!/usr/bin/python

from __future__ import division
from pylab import *
import sys
#if len(sys.argv) != 2:
#    print("Usage:  " + sys.argv[0] + " out-filename.pdf")
#    exit(1)

#pdffilename = sys.argv[1]

figure(1)
title('$g_{\sigma}$')
axvline(x=1, color='k', linestyle=':')
axhline(y=1, color='k', linestyle=':')
figure(2)
title('density')
figure(3)
title('nA')

def read_gs(ff):
  filename = "figs/wallsWB-0.%d0.dat" % (10*ff)
  print 'Using', filename
  data = loadtxt(filename)
  r = data[:,0]
  density = data[:,1]
  gsigma = data[:,2]
  nA = data[:,3]
  return r, density, gsigma, nA

ff = [.1, .2, .3, .4, .5]
for i in ff:
  r, density, gsigma, nA = read_gs(i)
  figure(1)
  plot(r, gsigma)
  figure(2)
  plot(r-3, density)
  figure(3)
  plot(r, nA)

show()
