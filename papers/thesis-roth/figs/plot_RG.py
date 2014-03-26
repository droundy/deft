from __future__ import division
import RG
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
amps = np.linspace(0,1,1000)
T = 0.1

def plotID_integrand(x):
    plt.plot(x,RG.integrand_ID(T,x,0))
    plt.ylabel('integrand of ID')
    plt.xlabel('amplitude of wave-packet')


def plotID_ref_integrand(x):
    plt.plot(x,RG.integrand_ID_ref(x))
    plt.ylabel('integrand of ID_ref')
    plt.xlabel('amplitude of wave-packet')

def plotphi(T,n,npart,i):
    plt.plot(n*RG.sigma**3*np.pi/6,RG.phi(T,n,npart,i))
    plt.ylabel(r'$\phi$')
    plt.xlabel(r'$\eta$')

# plotID_integrand(amps)
# plotID_ref_integrand(amps)

n = np.linspace(0,0.1,1000)

plotphi(T,n,0.171/(RG.sigma**3*np.pi/6),0)

plt.show()
