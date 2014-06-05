from __future__ import division
import numpy as np
import sys,matplotlib
import RG

# Make plot from Forte's data

forte_data = np.loadtxt('figs/forte_data.dat')

mrho = 0.8/(504-94)
brho = -mrho*94

mT = (0.65-1.85)/(345-52)
bT = 0.65-mT*345

rhostar = forte_data[:,0]*mrho + brho
T = forte_data[:,1]*mT + bT

eta = rhostar*np.pi/6

my_data0 = np.loadtxt('figs/npart_RG-i0-out.dat')
my_data1 = np.loadtxt('figs/npart_RG-i1-out.dat')

T0 = my_data0[:,0]
etavapor0 = my_data0[:,1]*np.pi*RG.sigma**3/6
etaliquid0 = my_data0[:,2]*np.pi*RG.sigma**3/6

T1 = my_data1[:,0]
etavapor1 = my_data1[:,1]*np.pi*RG.sigma**3/6
etaliquid1 = my_data1[:,2]*np.pi*RG.sigma**3/6

colors = np.array(['b-','g-','r-'])

if __name__ == '__main__':
  if 'show' not in sys.argv:
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  plt.plot(etavapor0, T0, colors[0],label='RG '+r'$i=0$')
  plt.plot(etaliquid0, T0, colors[0])

  plt.plot(etavapor1, T1, colors[1],label='RG '+r'$i=1$')
  plt.plot(etaliquid1, T1, colors[1])

  plt.plot(eta,T,colors[2],label='Forte 2011')

  plt.xlabel(r'$\eta$')
  plt.ylabel('T')
  plt.legend(loc=0)
  plt.title('l-v coexistence '+r'$\lambda_{SW}=1.5$')
#  plt.show()
  plt.savefig('figs/coexistence-RG-forte.pdf')
