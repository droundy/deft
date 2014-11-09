#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib, sys
if not "show" in sys.argv:
    matplotlib.use('Agg')

import pylab, numpy, os, glob
from pylab import pi

import styles

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + ' filling-fraction')
    exit(1)

rho = int(sys.argv[1])
#arg rho = [10, 20, 30, 40, 50, 60, 70, 80, 90]

pylab.figure()
data = []
names = []
lines = []
# eventually we will want to include this loadtx("figs/walls.dat") # just so things get rebuilt

for kT in [0.2, 1.0]:
    # input: ["figs/new-data/soft-wall-%04.2f-%04.2f.dat" % (rho*0.01, kT) for kT in [0.2, 1.0]]
    names.append('new kT = %g' % kT)
    fname = "figs/new-data/soft-wall-%04.2f-%04.2f.dat" % (rho*0.01, kT)
    data.append(numpy.loadtxt(fname))
    lines.append('-')

pylab.plot(data[0][:,0], data[0][:,2]*0.1, 'r:', label='$%g V_{wall}/kT$' % (0.1))

for i in range(len(data)):
    print 'plotting', names[i]
    pylab.plot(data[i][:,0], data[i][:,1], lines[i], label=names[i])

pylab.title('reduced density = %g' % (rho/100.0))
pylab.xlabel('$z/R$')
pylab.ylabel('reduced density')
pylab.legend(loc = 'best')

pylab.savefig('figs/soft-walls-%02d.pdf' % (rho))
pylab.show()
