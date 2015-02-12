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

def smooth(x, N):
    '''
    smooth(x,N) takes a 2D array x that has many columns, and averages
    out N nearby points.
    '''
    n0 = len(x)
    n = n0 - n0 % N
    x = x[:n]
    y = numpy.zeros_like(x[0::N])
    for i in range(N):
        y += x[i::N]
    return y/N

def average_positive_and_negative(data):
    for i in range(len(data)//2):
        data[i,1] = data[abs(data[:,0]) == abs(data[i,0]),1].mean()
    return data[:len(data)//2, :]
    

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
    lines.append('--')

pylab.plot(data[0][:,0], data[0][:,2]*0.1, 'r:', label='$V_{wall}$ (arbitrary units)')
for i in range(len(data)):
    print 'plotting', names[i]
    pylab.plot(data[i][:,0], data[i][:,1], lines[i], label=names[i])

N = 1
# softwall = numpy.loadtxt('figs/mc-soft-wall-0.5000-0.0100-2218.dat')
# pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N))

softwall = average_positive_and_negative(numpy.loadtxt('figs/mc-soft-wall-0.5000-0.2000-2195.dat'))
pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N), 'b-')
# softwall = average_positive_and_negative(numpy.loadtxt('figs/mc-soft-wall-0.5000-0.2000-2218.dat'))
# pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N), 'b-')
#softwall = numpy.loadtxt('figs/mc-soft-wall-0.5000-0.2000-2195.dat')
#pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N), 'b:')


# softwall = numpy.loadtxt('figs/mc-soft-wall-0.5000-0.5000-2218.dat')
# pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N))
softwall = average_positive_and_negative(numpy.loadtxt('figs/mc-soft-wall-0.5000-1.0000-2218.dat'))
pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N), 'g-')
# softwall = numpy.loadtxt('figs/mc-soft-wall-0.5000-2.5000-2218.dat')
# pylab.plot(smooth(15-abs(softwall[:,0]),N), smooth(softwall[:,1]/2**(-5.0/2.0),N))

pylab.title('reduced density = %g' % (rho/100.0))
pylab.xlabel('$z/R$')
pylab.ylabel('reduced density')
pylab.legend(loc = 'best')

pylab.xlim(0, 5)

pylab.savefig('figs/soft-walls-%02d.pdf' % (rho))
pylab.show()
