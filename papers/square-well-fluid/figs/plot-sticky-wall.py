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
lenx = float(sys.argv[3])
#arg lenx = [50]
lenyz = float(sys.argv[4])
#arg lenyz = [10]
T = float(sys.argv[5])
#arg T = [1.0]

plt.figure()

density, x = readandcompute.density_x('data/lv/ww%.2f-ff%.2f-%gx%g' % (ww,ff,lenyz,lenx), T)
plt.plot(x/2, density)

plt.xlabel(r'$z/\sigma$')
plt.ylabel(r'$\eta$')
plt.title(r'$\eta(z)$ with $\lambda = %g$, $\eta=%g$, and $T/\epsilon = %g$'
          % (ww, ff, T))

plt.savefig('figs/liquid-vapor-ww%.2f-ff%.2f-%gx%g-T%.2g.pdf' % (ww,ff,lenyz,lenx,T))

plt.show()
