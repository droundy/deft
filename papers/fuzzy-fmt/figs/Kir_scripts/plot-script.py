#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt

os.system('rm newmeltdataout.dat')

for gwidth in [.46, .459, .458, .4578, .4576, .4574, .4572, .457, .456]:
    os.system('../new-melting.mkdat 1.3 .1 %g 2' % (gwidth))

# this is a comment (not) describing the loadtxt.
data = np.loadtxt('newmeltdataout.dat')
plt.plot(data[:,0], data[:,1])

plt.title('Awesome melting plot!')
plt.xlabel('x label')

plt.show()
