from __future__ import division
import rg
import numpy as np
import matplotlib.pyplot as plt

# Plot f vs n
amps = np.linspace(0,1,1000)
T = 0.1

def plotID_integrand(x):
    plt.plot(x,rg.integrand_ID(x))
    plt.ylabel('integrand of ID')
    plt.xlabel('amplitude of wave-packet')
    plt.show()


def plotID_ref_integrand(x):
    plt.plot(x,rg.integrand_ID_ref(x))
    plt.ylabel('integrand of ID_ref')
    plt.xlabel('amplitude of wave-packet')
    plt.show()
