#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab
import SW
import RG

eta_conv = SW.sigma**3*pylab.pi/6

n = pylab.linspace(0,0.2,10000)
T = 0.3
npart = 0.0462316928818


pylab.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')
pylab.plot(n*eta_conv,RG.phi(T,n,npart,0),'r-',label='SW+RG; i=0')
pylab.plot(n*eta_conv,RG.phi(T,n,npart,1),'b--',label='SW+RG; i=1')

# pylab.plot(n*eta_conv,SW.f(T,n),'g-',label='SW')
# pylab.plot(n*eta_conv,RG.ftot(T,n,0),'r--',label='SW+RG; i=0')
# pylab.plot(n*eta_conv,RG.ftot(T,n,1),'b--',label='SW+RG; i=1')

pylab.title('T = %0.2f'%T)
pylab.ylabel(r'$\phi$')
# pylab.ylabel(r'$f$')
pylab.xlabel(r'$\eta$')
pylab.legend(loc=0)

pylab.savefig('figs/SW-RG-compare.pdf')
pylab.show()
