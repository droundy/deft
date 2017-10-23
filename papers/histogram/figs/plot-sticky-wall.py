#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readnew

if len(sys.argv) < 5:
    print("Usage: python {} 1.3 0.22 100 10".format(sys.argv[0]))
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]

plt.figure()

Ts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 10.0]

colors = { 0.1: 'r',
           0.5: 'k',
           0.6: 'y',
           0.7: 'm',
           0.8: 'g',
           0.9: 'b',
           1.0: 'r',
           5.0: 'k',
           10.0: 'c',
       }
def color(T):
    try:
        return colors[T]
    except:
        return ''
lines = ['-', '--', ':','-.', '.']

first_method = True
the_first_method = ''

methods = [ '-tmi3', '-toe3', '-tmmc', '-wltmmc-0.8-1e-10'] #, '-tmi']
first_temperature = [True]*len(methods)

for i in range(len(methods)):
    method = methods[i]
    fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s' % (ww,ff,lenx,lenyz,method)
    fname = fbase + '-density.dat'
    try:
        minT = readnew.minT(fname)
        convergedT = readnew.convergedT(fname)

        for T in Ts:
            if T >= minT and T >= convergedT*1.0:
                density, x = readnew.density_x(fbase, T)
                plt.plot(x/2, density, color(T)+lines[i])
                if first_method or method == the_first_method:
                    if first_temperature[i]:
                        plt.plot(x/2, density, color(T)+lines[i],
                                 label='T=%g %s (converged to %.2g)' % (T, method[1:], convergedT))
                        first_temperature[i] = False
                    else:
                        plt.plot(x/2, density, color(T)+lines[i], label='T=%g' % T)
                    the_first_method = method
                    first_method = False
                elif first_temperature[i]:
                    plt.plot(x/2, density, color(T)+lines[i],
                             label='T=%g %s (converged to %.2g)' % (T, method[1:], convergedT))
                    first_temperature[i] = False
                else:
                    plt.plot(x/2, density, color(T)+lines[i])
    except:
        pass

plt.ylim(0)
plt.xlabel(r'$z/\sigma$')
plt.ylabel(r'$\eta$')
plt.legend(loc='best')
plt.title(r'$\eta(z)$ with $\lambda = %g$ and $\eta=%g$' % (ww, ff))

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g.pdf' % (ww,ff,lenx,lenyz))

plt.show()
