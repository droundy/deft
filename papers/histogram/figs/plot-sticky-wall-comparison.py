#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readnew
import glob
import os


if len(sys.argv) < 5:
    print("Usage: python {} 1.3 0.22 100 0.8".format(sys.argv[0]))
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]
T = float(sys.argv[5])
#arg T = [0.8]

plt.figure()

methods = [ 'tmi3', 'tmmc', '-toe3', '-wltmmc-1-0.0001']

fname = glob.glob('data/lv/ww%.2f-ff%.2f-%gx%g*-density.dat' % (ww,ff,lenx,lenyz))
fbase = sorted(fname, key=os.path.basename)

for i in range(len(fbase)):
    density, x = readnew.density_x(fbase[i], T)

    plt.plot(x/2, density,label='T=%g-%gx%g%s' % (T, lenx,lenyz, methods[i]))

plt.ylim(0)
plt.xlabel(r'$z/\sigma$')
plt.ylabel(r'$\eta$')
plt.legend(loc='best')
plt.title(r'$\eta(z)$ with $\lambda = %g$ and $\eta=%g$' % (ww, ff))

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g-%g.pdf' % (ww,ff,lenx,lenyz,T))

plt.show()
