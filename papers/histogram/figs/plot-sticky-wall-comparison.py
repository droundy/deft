#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readnew

if len(sys.argv) < 8:
    print("Usage: python {} 1.3 0.22 100 10 0.8".format(sys.argv[0]))
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]
lenyz_ref1 = float(sys.argv[5])
#arg lenyz = [10]
lenyz__ref2 = float(sys.argv[6])
#arg lenyz = [10]
T = float(sys.argv[7])
#arg T = [0.8]

plt.figure()

methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']
method = methods[0]

fbase_1 = 'data/lv/ww%.2f-ff%.2f-%gx%g%s' % (ww,ff,lenx,lenyz_ref1,method)
fbase_2 = 'data/lv/ww%.2f-ff%.2f-%gx%g%s' % (ww,ff,lenx,lenyz__ref2,method)
    
fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s' % (ww,ff,lenx,lenyz,method)

fname_1 = fbase_1 + '-density.dat'
fname_2 = fbase_2 + '-density.dat'
fname = fbase + '-density.dat'

density_1, x_1 = readnew.density_x(fbase_1, T)
density_2, x_2 = readnew.density_x(fbase_2, T)
density, x = readnew.density_x(fbase, T)

plt.plot(x_1/2, density_1,label='T=%g-%gx%g%s' % (T, lenx, lenyz_ref1, method[0:]))
plt.plot(x_2/2, density_2,label='T=%g-%gx%g%s' % (T, lenx, lenyz__ref2, method[0:]))
plt.plot(x/2, density,label='T=%g-%gx%g%s' % (T, lenx, lenyz, method[0:]))

plt.ylim(0)
plt.xlabel(r'$z/\sigma$')
plt.ylabel(r'$\eta$')
plt.legend(loc='best')
plt.title(r'$\eta(z)$ with $\lambda = %g$ and $\eta=%g$' % (ww, ff))

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g-%g.pdf' % (ww,ff,lenx,lenyz,T))

plt.show()
