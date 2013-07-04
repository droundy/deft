#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
#matplotlib.use('Agg')
import pylab, numpy, sys
import os.path
import matplotlib.colors as mcolors

ff = .3
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])

def read_walls(ff):
    filename = "figs/mc/wallsMC-pair-%02.1f-density.dat" % ff
    if (os.path.isfile(filename) == False):
        print "File %s does not exist" %filename
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

def read_gs(ff):
  filename = "figs/wallsWB-%03.2f.dat" % ff
  if (os.path.isfile(filename) == False):
    print "File %s does not exist" %filename
    sys.exit(1)
  data = pylab.loadtxt(filename)
  r = data[:,0]
  density = data[:,1]
  gsigma = data[:,2]
  nA = data[:,3]
  return r, density

den = read_walls(ff)
r, den2 = read_gs(ff)


N = sum(den[:,1])/len(den[:,1])*30*30*30
Ndft = sum(den2)/(len(den2)-3/.01)*30*30*30
print N, Ndft


pylab.xlim(0,15)
pylab.plot(den[:,0], den[:,1])
pylab.plot(r-3, den2)
pylab.title('density(z), $ff = %g$' %(ff))
pylab.xlabel("z")
pylab.ylabel("density (balls/unit$^3$)")
pylab.legend(['mc', 'deft'], loc='best')
pylab.show()
