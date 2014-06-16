from __future__ import division
import numpy as np
import sys,matplotlib
import RG

# Make plot from Forte's data

forte_data = np.loadtxt('figs/forte_data.dat') # not converged
forte_data_conv = np.loadtxt('figs/forte_data_conv.dat') # Converged
forte_data_MD = np.loadtxt('figs/forte_data_MD.dat') # MD

mrho = 0.8/(504-94)
brho = -mrho*94

mT = (0.65-1.85)/(345-52)
bT = 0.65-mT*345

### Not converged ###
rhostar_forte = forte_data[:,0]*mrho + brho
T_forte = forte_data[:,1]*mT + bT + 0.21
eta_forte = rhostar_forte*np.pi/6

### Converged ###
rhostar_conv = forte_data_conv[:,0]*mrho + brho
T_forte_conv = forte_data_conv[:,1]*mT + bT
eta_conv = rhostar_conv*np.pi/6

### MD sim ###
rhostar_MD = forte_data_MD[:,0]*mrho + brho
T_forte_MD = forte_data_MD[:,1]*mT + bT + 0.22
eta_MD = rhostar_MD*np.pi/6

### My RGT ###
my_data0 = np.loadtxt('figs/npart_RG-i0-out.dat')
# my_data1 = np.loadtxt('figs/npart_RG-i1-out.dat')

T0 = my_data0[:,0]
etavapor0 = my_data0[:,1]*np.pi*RG.sigma**3/6
etaliquid0 = my_data0[:,2]*np.pi*RG.sigma**3/6

# T1 = my_data1[:,0]
# etavapor1 = my_data1[:,1]*np.pi*RG.sigma**3/6
# etaliquid1 = my_data1[:,2]*np.pi*RG.sigma**3/6

colors = np.array(['b-','g-','ro','c-'])

if __name__ == '__main__':
  if 'show' not in sys.argv:
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  plt.plot(etavapor0, T0, colors[0],label='RG '+r'$i=0$')
  plt.plot(etaliquid0, T0, colors[0])

  # plt.plot(etavapor1, T1, colors[1],label='RG '+r'$i=1$')
  # plt.plot(etaliquid1, T1, colors[1])

  plt.plot(eta_forte,T_forte,colors[1],label='Forte 2011')

  plt.plot(eta_MD,T_forte_MD,colors[2],label='Elliot 1999')

  plt.plot(eta_conv,T_forte_conv,colors[3],label='Forte 2011, converged RGT')

  plt.xlabel(r'$n$')
  plt.ylabel(r'$T$')
  plt.legend(loc=0)
  plt.title('l-v coexistence '+r'$\lambda_{SW}=1.5$')
#  plt.show()
  plt.savefig('figs/coexistence-RG-forte.pdf')
