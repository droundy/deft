#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, scipy

from scipy.interpolate import interp1d

if len(sys.argv) != 4:
    print("Usage:  " + sys.argv[0] + " cell-width radius filling-fraction")
    exit(1)

width = float(sys.argv[1])
radius = float(sys.argv[2])
ff = float(sys.argv[3])

if width == 24 and radius == 2:
    num = [0,
           332,
           501,
           659,
           751,
           987,
           1153,
           1320,
           1647]
    ffvals = [0,
              0.1008,
              0.1515,
              0.2,
              0.2280,
              0.2995,
              0.350,
              0.4,
              0.5]
elif width == 32 and radius == 6:
    num = [0,
           764,
           1410,
           1528,
           2291,
           3055]
    ffvals = [0,
              0.1,
              0.1846,
              0.2,
              0.299,
              0.398]
else:
    print "I don't have enough data"
    sys.exit(1)

# interpolate beyond highest...
num.append(num[-1]*100)
ffvals.append(ffvals[-1]*100)
num = numpy.array(num)
ffvals = numpy.array(ffvals)

print numpy.interp(ff, ffvals, num)

print interp1d(ffvals, num, kind='cubic')(ff)
