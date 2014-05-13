#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import RG

# # define the colors/symbols
# colors = np.array(['b-','g-','r-'])

# # normal
# plt.figure()
# for i in range(1):
#   # Read in data
#   data = np.loadtxt('figs/npart_RG-i%d-out.dat'%i)

#   T = data[:,0]
#   etavapor = data[:,1]*np.pi*RG.sigma**3/6
#   etaliquid = data[:,2]*np.pi*RG.sigma**3/6

#   # Plot the curve
#   plt.plot(etavapor, T, colors[i],label='RG '+r'$i=$ '+'%d'%i)
#   plt.plot(etaliquid, T, colors[i])

# plt.xlabel(r'$\eta$')
# plt.ylabel('T')
# plt.legend(loc=0)
# plt.title('Liquid-Vapor Coexistence')
