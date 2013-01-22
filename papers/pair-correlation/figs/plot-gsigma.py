#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
#matplotlib.use('Agg')

import pylab, numpy, sys
#if len(sys.argv) != 2:
#    print("Usage:  " + sys.argv[0] + " out-filename.pdf")
#    exit(1)

#pdffilename = sys.argv[1]

pylab.figure(1)
pylab.title('$g_{\sigma}$')
pylab.xlabel('r')
pylab.axvline(x=1, color='k', linestyle=':')
pylab.axhline(y=1, color='k', linestyle=':')

def read_gs(ff):
    filename = "wallsWB-0.%d0.dat" % (10*ff)
    print 'Using', filename
    data = numpy.loadtxt(filename)
    r = data[:,0]
    gsigma = data[:,2]
    return r, gsigma

ff = [.1, .2, .3, .4, .5]
for i in ff:
    r, gsigma = read_gs(i)

    pylab.plot(r, gsigma)


pylab.show()
