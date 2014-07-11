from __future__ import division
import RG_v2 as RG
import RG as RG1
import numpy as np
import matplotlib.pyplot as plt

# Conversions and limits
eta_conv = RG.sigma**3*np.pi/6
nmax = 0.75/eta_conv
ns = nmax*np.exp(np.arange(-15, 0, 1e-3))

# Plots
def plotphi(T,n,mu,i,label):
  plt.plot(n*eta_conv,RG.phi(T,n,mu,i),label=label)

def plotphi_old(T,n,npart,i,label):
  plt.plot(n*eta_conv,RG1.phi(T,n,npart,i),label=label)

def plotPress(T,n,i):
  plt.plot(n*eta_conv,RG.P(T,n,i))

### Plot stuff ###
T = 0.5

mu = -7
plotphi(T,ns,mu,1,label='mu=%.3f'%mu)

mu = -6
plotphi(T,ns,mu,1,label='mu=%.3f'%mu)

mu = -5
plotphi(T,ns,mu,1,label='mu=%.3f'%mu)

mu = -4
plotphi(T,ns,mu,1,label='mu=%.3f'%mu)

### Titles, etc. ###
plt.title('Grand free energy per unit volume; T=%.2f'%T)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\phi(T,n,\mu)$'+'/SW units')

# plt.xlim(0.0,0.3)
# plt.ylim(-0.02,0.0)
plt.legend(loc=0)

plt.savefig('figs/mu_comparison-i1.pdf')
#plt.show()
