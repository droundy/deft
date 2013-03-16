#!/usr/bin/python

from __future__ import division

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

from pylab import *

dx = 0.01
x = arange(0, 1 + dx/2, dx)
harmonic = (1-x)**2

figure(figsize=(4.5,4))
plot(x,harmonic)

xlabel('$r/R$')
ylabel('$V/V_{max}$')

savefig('figs/harmonic.pdf')
