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

n_reduced = float(sys.argv[1])
#arg n_reduced = [0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5, 2.0]

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
        y += x[i::N, :]
    return y/N

pylab.figure()
data = []
names = []
lines = []
# eventually we will want to include this loadtx("figs/walls.dat") # just so things get rebuilt

for kT in [0.2, 0.5, 1.0]:
    # input: ["figs/mcwalls-%.4f-%.4f-*.dat" % (n_reduced, kT) for kT in [0.2, 0.5, 1.0]]
    print 'looking for', 'figs/mcwalls-%.4f-%.4f*.dat' % (n_reduced, kT)
    for fname in glob.glob('figs/mcwalls-%.4f-%.4f*.dat' % (n_reduced, kT)):
        print 'examining', fname
        d = numpy.loadtxt(fname)
        d = smooth(d, 10)
        pylab.plot(15-abs(d[:,0]), d[:,1]/2.0**(-5.0/2.0), styles.mc[kT], label=fname)

# names.append('DFT kT = 0')
# data.append(numpy.loadtxt("figs/wallshard-%.4f-%.2f.dat" % (0.0, n_reduced)))
# lines.append(styles.dft[0])

for kT in [0.2, 0.5, 1.0]:
    # input: ["figs/new-data/wall-%04.2f-%04.2f.dat" % (n_reduced, kT) for kT in [0.2, 0.5, 1.0]]
    names.append('new kT = %g' % kT)
    fname = "figs/new-data/wall-%04.2f-%04.2f.dat" % (n_reduced, kT)
    data.append(numpy.loadtxt(fname))
    lines.append(styles.new_dft_code[kT])

# for kT in [0.001, 0.01, 0.03]:
#     # input: "figs/wallssoft-0.0010-%.2f.dat" % (n_reduced)
#     # input: "figs/wallssoft-0.0100-%.2f.dat" % (n_reduced)
#     # input: "figs/wallssoft-0.0300-%.2f.dat" % (n_reduced)
#     fname = "figs/wallssoft-%.4f-%.2f.dat" % (kT, n_reduced)
#     if os.path.exists(fname):
#         names.append('DFT kT = %4.2f' % kT)
#         data.append(numpy.loadtxt(fname))
#         lines.append(styles.dft[kT])
#         print 'found', fname
#     else:
#         print fname, 'does not exist'


for i in range(len(data)):
    print 'plotting', names[i]
    pylab.plot(data[i][:,0], data[i][:,1], lines[i], label=names[i])

pylab.title('Reduced Density = %f' % (n_reduced))
pylab.xlabel('$z/R$')
pylab.ylabel('Reduced density')
pylab.legend(loc = 'best')
#pylab.xlim(-0.1, 15)
#pylab.ylim(.25, .6)
pylab.savefig('figs/walls-%02.0f.pdf' % (n_reduced*100))
pylab.show()
