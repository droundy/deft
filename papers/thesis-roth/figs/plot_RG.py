from __future__ import division
import RG
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
R = 1 # HS radius
eta_conv = RG.sigma**3*np.pi/6
nmax = 0.6/eta_conv
n = nmax*np.exp(np.arange(-15, 0, 1e-3))

def plotID_integrand(x):
  plt.plot(x,RG.integrand_ID(T,x,0))
  plt.ylabel('integrand of ID')
  plt.xlabel('amplitude of wave-packet')


def plotID_ref_integrand(x):
  plt.plot(x,RG.integrand_ID_ref(x))
  plt.ylabel('integrand of ID_ref')
  plt.xlabel('amplitude of wave-packet')

def plotphi(T,n,npart,i):
  plt.plot(n*eta_conv,RG.phi(T,n,npart,i))
  plt.ylabel(r'$\phi$')
  plt.xlabel(r'$\eta$')

def plotPress(T,n,i):
  plt.plot(n*eta_conv,RG.P(T,n,i))
  plt.ylabel('P')
  plt.xlabel(r'$\eta$')

# plotphi(0.1,n,0.025,0)

plotPress(0.1,n,0)
plotPress(0.1,n,1)

plt.show()
