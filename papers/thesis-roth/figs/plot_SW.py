from __future__ import division
import SW
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
R = 1 # HS radius
#nmax = np.sqrt(2)/R # max density is that of close-packed spheres of radius R
nmax = 0.4/(SW.sigma**3*np.pi/6)
n = nmax*np.exp(np.arange(-15, 0, 1e-3))

def plotPhi(T,npart):
    plt.plot(n*SW.sigma**3*np.pi/6,SW.phi(T,n,npart))
    plt.ylabel(r'$\phi$ (SW units)')

    plt.xlabel(r'$\eta$')

def plotF(T):
    plt.subplot(211)
    plt.plot(n,SW.f(T,n))
    plt.ylabel(r'$f$ (SW units)')
    plt.title('T=%0.2f'%T)

    plt.subplot(212)
    plt.plot(n,SW.df_dn(T,n))
    plt.ylabel(r'$\frac{df}{dn}$ (SW units)')

#    plt.savefig('f_T%0.2.pdf'%T)

# plt.show()

T = 0.01
plotPhi(T, 0.05/(SW.sigma**3*np.pi/6))
plt.savefig('SW-T%0.2f.pdf'%T)
plt.show()
