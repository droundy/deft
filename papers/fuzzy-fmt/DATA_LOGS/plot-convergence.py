#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('DATA_LOGS/n1.3_gw_0.01.dat')

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend(loc='best')

plt.show()
