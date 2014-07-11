from __future__ import division

import matplotlib,sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import RG_v2 as RG

def correction(Temp,num,chempot,i1,i2):
    return RG.phi(Temp,num,chempot,i1) - RG.phi(Temp,num,chempot,i2)


ns = np.arange(1e-10,0.6/(np.pi*RG.sigma**3/6),1e-3/(np.pi*RG.sigma**3/6))

T = 1.31705231562
mu = -6.21486388264

T = 1.18142906562
mu = -6.02699427421

T = 1.2219968
mu = -6.08764244458

T = 0.5
mu = -5.75178857763

#plt.plot(ns*np.pi*RG.sigma**3/6, correction(T,ns,mu,1,0), label=r'$(T,\mu)=$ (%g,%g) difference'%(T,mu))
plt.plot(ns*np.pi*RG.sigma**3/6, RG.phi(T,ns,mu,0), label=r'$\phi(%g,\eta,%g,0)$'%(T,mu))
plt.plot(ns*np.pi*RG.sigma**3/6, RG.phi(T,ns,mu,1), label=r'$\phi(%g,\eta,%g,1)$'%(T,mu))
plt.plot(ns*np.pi*RG.sigma**3/6, RG.phi(T,ns,mu,2), label=r'$\phi(%g,\eta,%g,2)$'%(T,mu))

plt.xlabel(r'$\eta$')
plt.ylabel(r'$\phi(T,\eta,\mu,i)$')
plt.legend(loc=0)

plt.savefig('figs/RG_comparison.pdf')
# plt.show()
