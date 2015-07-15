#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readandcompute

ww = float(sys.argv[1])
#arg ww = [1.3]
kTslices = eval(sys.argv[2])
#arg kTslices = [[1,2,3,4,5,6]]

plt.figure()

data = np.loadtxt('data/coexistence/ww%g.dat' % ww)
plt.plot(data[:,0], data[:,1]*np.pi*4/3)
plt.plot(data[:,0], data[:,2]*np.pi*4/3)
plt.plot(data[:,0], data[:,3]*np.pi*4/3)

#input: ['data/coexistence/ww%g-kT%g.dat' % (ww, kT) for kT in kTslices]
for kT in kTslices:
    fname = 'data/coexistence/ww%g-kT%g.dat' % (ww, kT)
    ugh = np.loadtxt(fname)
    plt.plot(kT + ugh[:,1] - ugh[:,1].min(), ugh[:,0]*np.pi*4/3, ':')

plt.xlabel(r'$T$')
plt.ylabel(r'$\eta$')

plt.savefig('figs/coexistence-ww%s.pdf' % (('%g'%ww).replace('.','_')) )

plt.show()
