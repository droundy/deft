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
T = float(sys.argv[4])
#arg T = [0.8]

plt.figure()

methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']
method = methods[0]
 
fname = glob.glob('data/lv/ww%.2f-ff%.2f-%gx*%s-density.dat' % (ww,ff,lenx,method))
fbase = sorted(fname, key=os.path.getsize)

for i in range(len(fbase)):
    density, x = readnew.density_x(fbase[i], T)
    
    # this seems like a clunkey way to correctly label the file
    tag = fbase[i].split("x",1)[1]
    lenyz = float(tag.split("-",1)[0])
    
    plt.plot(x/2, density,label='T=%g-%gx%g%s' % (T, lenx,lenyz, method[0:]))

plt.ylim(0)
plt.xlabel(r'$z/\sigma$')
plt.ylabel(r'$\eta$')
plt.legend(loc='best')
plt.title(r'$\eta(z)$ with $\lambda = %g$ and $\eta=%g$' % (ww, ff))

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%g-%g.pdf' % (ww,ff,lenx,T))

plt.show()
