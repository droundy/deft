#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib, sys
if not "show" in sys.argv:
    matplotlib.use('Agg')

import pylab, numpy, os, glob
from pylab import pi

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + ' filling-fraction')
    exit(1)

ff = int(sys.argv[1])
#arg ff = [10, 20, 30, 40]

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

pylab.figure()
data = []
names = []
# eventually we will want to include this loadtx("figs/walls.dat") # just so things get rebuilt

for kT in [0.0, 0.1, 0.01, 0.001, 0.0001]:
    # input: "figs/mcwalls-0.%02d00-*.dat" % (ff)
    print 'looking for', 'figs/mcwalls-0.%02d00-%.4f*.dat' % (ff, kT)
    for fname in glob.glob('figs/mcwalls-0.%02d00-%.4f*.dat' % (ff, kT)):
        print 'examining', fname
        d = numpy.loadtxt(fname)
        d = smooth(d, 1)
        pylab.plot(abs(d[:,0]), d[:,1]*(4*pi/3), label=fname)

names.append('DFT kT = 0')
data.append(numpy.loadtxt("figs/wallshard-%.4f-%.2f.dat" % (0.0, ff*0.01)))

for kT in [0.01, 0.02, 0.03]:
    # eventually: nput: "figs/walls-0.%02d-T*.dat" % (ff)
    fname = "figs/wallssoft-%.4f-%.2f.dat" % (kT, ff*0.01)
    if os.path.exists(fname):
        names.append('DFT kT = %4.2f' % kT)
        data.append(numpy.loadtxt(fname))
        print 'found', fname
    else:
        print fname, 'does not exist'
for i in range(len(data)):
    print 'plotting', names[i]
    pylab.plot(data[i][:,0]-1.5, data[i][:,1], label=names[i])

pylab.title('Packing fraction = %f' % (ff/100.0))
pylab.xlabel('$z/R$')
pylab.ylabel('Local Filling Fraction')
pylab.legend(loc = 'best')
pylab.xlim(-0.1, 5)
#pylab.ylim(.25, .6)
pylab.savefig('figs/walls-%02d.pdf' % (ff))
pylab.show()
