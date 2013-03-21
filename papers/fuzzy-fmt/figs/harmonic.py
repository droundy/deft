#!/usr/bin/python

from __future__ import division

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

from pylab import *

dx = 0.01
x = arange(0, 2.1 + dx/2, dx)
harmonic = (2-x)**2/4
harmonic[x>2] = 0

figure(figsize=(4.5,4))
plot(x,harmonic)
xlim(0,2.2)

xlabel('$r/R$')
ylabel('$V/V_{max}$')

savefig('figs/harmonic.pdf')
