#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab as plt
import SW
import RG

eta_conv = SW.sigma**3*plt.pi/6

n = plt.linspace(1e-8,0.2,10000)

### Low temp ###

T = 0.8
npart = 0.04

plt.figure()
plt.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
plt.ylim(-0.01,0.04)
plt.xlim(0,0.45)
#plt.legend(loc=0)

plt.savefig('figs/SW-phi-lowT.pdf')


### High Temp ###

T = 1.32
npart = 0.0360850913603

plt.figure()
plt.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
#plt.legend(loc=0)

plt.savefig('figs/SW-phi-highT.pdf')
