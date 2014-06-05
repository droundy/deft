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

T = 0.3
npart0 = 0.0462316928818
npart1 = 0.0776670124628

plt.figure()
plt.plot(n*eta_conv,RG.phi(T,n,npart0,0),label='i=0')
plt.plot(n*eta_conv,RG.phi(T,n,npart1,1),label='i=1')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
plt.ylim(-0.04,0.1)
plt.xlim(0,0.7)
plt.legend(loc=0)

plt.savefig('figs/phi-RG-lowT.pdf')

### High Temp ###

T = 1.32
npart0 = 0.0360850913603
npart1 = 0.116216647163

plt.figure()
plt.plot(n*eta_conv,RG.phi(T,n,npart0,0),label='i=0')
plt.plot(n*eta_conv,RG.phi(T,n,npart1,1),label='i=1')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
plt.ylim(-2,2)
plt.xlim(0,0.7)
plt.legend(loc=0)

plt.savefig('figs/phi-RG-highT.pdf')
