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

first_method = True
the_first_method = ''

methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']
first_temperature = [True]*len(methods)
method = methods[0]

#print('ww%.2f-ff%.2f-%gx%g%s-movie') % (ww,ff,lenx,lenyz,method)
fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s-movie/000070' % (ww,ff,lenx,lenyz,method)


lndos,ps = readnew.e_lndos_ps(fbase)

plt.semilogx(ps,lndos)

#plt.ylim(0)
#plt.xlabel(r'$Pessimistic Samples$')
#plt.ylabel(r'$lndos$')
#plt.legend(loc='best')
#plt.title(r'$lndos$ with $\lambda = %g$ and $\eta=%g$' % (ww, ff))

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-Error.pdf' % (ww,ff,lenx,lenyz,method))
plt.show()
