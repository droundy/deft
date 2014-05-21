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
npart = 0.0462316928818

plt.figure()
# plt.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')
plt.plot(n*eta_conv,RG.phi(T,n,npart,0),label='i=0')
plt.plot(n*eta_conv,RG.phi(T,n,npart,1),label='i=1')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
plt.legend(loc=0)

plt.savefig('figs/phi-RG-lowT.pdf')

# plt.figure()
# # plt.plot(n*eta_conv,SW.f(T,n),'g-',label='SW')
# plt.plot(n*eta_conv,RG.ftot(T,n,0),'r--',label='i=0')
# #plt.plot(n*eta_conv,RG.ftot(T,n,1),'b--',label='i=1')

# plt.title('T = %0.2f'%T)
# plt.ylabel(r'$f$')
# plt.xlabel(r'$\eta$')
# plt.legend(loc=0)

# I introduced a misspelling below, so scons won't think this should be built
# plt.save fig('figs/SW-RG-compare-f-lowT.pdf')
# plt.show()

### High Temp ###

T = 1.32
npart = 0.0360850913603

plt.figure()
# plt.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')
plt.plot(n*eta_conv,RG.phi(T,n,npart,0),label='i=0')
plt.plot(n*eta_conv,RG.phi(T,n,npart,1),label='i=1')

plt.title('T = %0.2f'%T)
plt.ylabel(r'$\phi$')
plt.xlabel(r'$\eta$')
plt.legend(loc=0)

plt.savefig('figs/phi-RG-highT.pdf')

# plt.figure()
# # plt.plot(n*eta_conv,SW.f(T,n),'g-',label='SW')
# plt.plot(n*eta_conv,RG.ftot(T,n,0),'r--',label='i=0')
# #plt.plot(n*eta_conv,RG.ftot(T,n,1),'b--',label='i=1')

# plt.title('T = %0.2f'%T)
# plt.ylabel(r'$f$')
# plt.xlabel(r'$\eta$')
# plt.legend(loc=0)

# I introduced a misspelling below, so scons won't think this should be built
# plt.save fig('figs/SW-RG-compare-f-highT.pdf')




# plt.show()
