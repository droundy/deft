from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import Hughes as H

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Oct 2013

def plotphi(T,prefac):

    # Plot P vs n
    nmax = 1.1*H.nw
    n = nmax*np.exp(np.arange(-15, 0, 1e-3))
    temp = str(T)

    # Plot P vs n
    # plt.title('T = ',temp)
    # plt.subplot(211)
    # plt.plot(n*H.conv_n,H.p(T,n)*H.conv_p)
    # plt.ylabel('P (atm)')

    # Plot GFE vs n
    # plt.subplot(212)
    plt.plot(n*H.conv_n,H.Phi(T,n,prefac)*H.conv_p)
    plt.ylabel(r'$\phi$ (atm)')

    plt.xlabel(r'$\rho$ (g/mL)')


def plotphi_alt(T,nparticular):
    nmax = 1.1*H.nw
    n = nmax*np.exp(np.arange(-15, 0, 1e-3))
    temp = str(T)

    plt.plot(n*H.conv_n,H.Phi_alt(T,n,nparticular)*H.conv_p)
    plt.ylabel(r'$\phi$ (atm)')

    plt.xlabel(r'$\rho$ (g/mL)')
