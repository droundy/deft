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
Ns = eval(sys.argv[3])
#arg Ns = [[500, 1372, 400]]
T = float(sys.argv[4])
#arg T = [1.0]

plt.figure()

for N in Ns:
    g, r = readandcompute.g_r('data/mc/ww%.2f-ff%.2f-N%d' % (ww,ff,N), T)
    plt.plot(r/2, g, '-', label='$N=%d$' % N)

plt.legend(loc='best')
plt.xlabel(r'$r/\sigma$')
plt.ylabel(r'$g(r)$ not properly normalized')
plt.title(r'$g(r)$ with $\lambda = %g$, $\eta=%g$, and $T/\epsilon = %g$'
          % (ww, ff, T))

plt.savefig('figs/radial-distribution-ww%.2f-ff%.2f-T%.2g.pdf' % (ww,ff,T))

plt.show()
