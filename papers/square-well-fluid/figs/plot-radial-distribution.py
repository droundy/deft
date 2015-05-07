#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readandcompute

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.3]
N = int(sys.argv[3])
#arg N = [30]
T = float(sys.argv[4])
#arg T = [1.0]

plt.figure()

g, r = readandcompute.g_r('data/mc/ww%.2f-ff%.2f-N%d' % (ww,ff,N), ff, T)

plt.plot(r/2, g)
plt.xlabel(r'$r/\sigma$')
plt.ylabel(r'$g(r)$ not properly normalized')
plt.title(r'$g(r)$ with $\lambda = %g$, $\eta=%g$, $N = %d$, and $T/\epsilon = %g$'
          % (ww, ff, N, T))

plt.savefig('radial-distribution-ww%.2f-ff%.2f-N%d-T%.2g.pdf' % (ww,ff,N,T))

plt.show()
