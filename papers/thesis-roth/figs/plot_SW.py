from __future__ import division
import SW
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
R = 1 # HS radius
nmax = np.sqrt(2)/R # max density is that of close-packed spheres of radius R
n = nmax*np.exp(np.arange(-15, 0, 1e-3))

def plotPhi(T,npart):
    plt.plot(n,SW.phi(T,n,npart))
    plt.ylabel(r'$\phi$ (SW units)')

    plt.xlabel('n (SW units)')

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

T = 3.85
plotPhi(T, 0.0670937256798)
plt.savefig('SW-T%0.2f.pdf'%T)
