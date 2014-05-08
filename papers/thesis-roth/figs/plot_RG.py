from __future__ import division
import RG
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
R = 1 # HS radius
eta_conv = RG.sigma**3*np.pi/6
nmax = 0.4/eta_conv
n = nmax*np.exp(np.arange(-15, 0, 1e-3))

def plotphi(T,n,npart,i):
  plt.plot(n*eta_conv,RG.phi(T,n,npart,i),label=r'$i=$'+'%d'%i)

def plotPress(T,n,i):
  plt.plot(n*eta_conv,RG.P(T,n,i))

plotphi(1.3299923125,n,0.0427515212938,0)
#plotphi(1.3299923125,n,0.0427515212938,1)

plt.title('Grand free energy per unit volume')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\phi(T,n,\mu)$'+'/SW units')

plt.xlim(0.0,0.3)
plt.ylim(-0.02,0.0)
plt.legend(loc=0)

plt.savefig('figs/plot-phi-RG.pdf')
