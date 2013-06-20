#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, os, glob
from pylab import pi

if len(sys.argv) != 1:
    print("Usage:  " + sys.argv[0])
    exit(1)

def smooth(x, N):
    '''
    smooth(x,N) takes a 2D array x that has many columns, and averages
    out N nearby points.
    '''
    n0 = len(x[:,0])
    n = n0 - n0 % N
    x = x[:n,:]
    y = numpy.zeros_like(x[0::N,:])
    for i in range(N):
        print 'y shape', y.shape
        print 'x[i::N,:] shape', x[i::N, :].shape
        y += x[i::N, :]
    return y/N

def plotff(ff):
    pylab.figure()
    data = []
    names = []
    for kT in [0.0, 0.01, 0.02, 0.03]:
        fname = "figs/walls-%04.2f-T%05.3f.dat" % (ff, kT)
        if os.path.exists(fname):
            names.append('kT = %4.2f' % kT)
            data.append(numpy.loadtxt(fname))
    for i in range(len(data)):
        pylab.plot(data[i][:,0]-3, data[i][:,1]*(4*pi/3), label=names[i])

    for kT in [0.0, 0.1, 0.01, 0.001, 0.0001]:
        print 'looking for', 'figs/mcwalls-%.4f-%.4f*.dat' % (ff, kT)
        for fname in glob.glob('figs/mcwalls-%.4f-%.4f*.dat' % (ff, kT)):
            print 'examining', fname
            d = numpy.loadtxt(fname)
            d = smooth(d, 10)
            pylab.plot(d[:,0], d[:,1]*(4*pi/3), label=fname)
    pylab.title('Packing fraction = %.1f' % ff)
    pylab.legend()
    pylab.xlim(xmax=14)
    pylab.savefig('figs/walls-%02d.pdf' % (ff*100))

for ff in [0.1, 0.2, 0.3, 0.4, 0.5]:
    plotff(ff)
