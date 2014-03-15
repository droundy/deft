from __future__ import division
import SW
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
R = 1 # HS radius
nmax = np.sqrt(2)/R # max density is that of close-packed spheres of radius R
n = nmax*np.exp(np.arange(-15, 0, 1e-3))
npart = 1e-1
T = 276

def plotPhi(T,npart):
    plt.plot(n,SW.phi(T,n,npart))
    plt.ylabel(r'$\phi$ (SW units)')

    plt.xlabel('n (SW units)')

#    plt.savefig(str(npart)+'b.png')

def plotF(T):
    plt.subplot(211)
    plt.plot(n,SW.f(T,n))
    plt.ylabel(r'$f$ (SW units)')
    plt.title(r'$T=$'+str(T))

    plt.subplot(212)
    plt.plot(n,SW.df_dn(T,n))
    plt.ylabel(r'$\frac{df}{dn}$ (SW units)')

#    plt.savefig('f_'+str(T)+'.png')

# plt.show()
