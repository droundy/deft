#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab as plt
import SW
import RG

eta_conv = SW.sigma**3*plt.pi/6

def plotphi(T, npart, filename, nv=0, nl=0, phi_min=0, phi_max=0):
  n = plt.linspace(1e-8,0.2,20000)
  plt.figure()
  plt.plot(n*eta_conv,SW.phi(T,n,npart),'g-',label='SW')

  plt.title('T = %0.2f'%T)
  plt.ylabel(r'$\phi$')
  plt.xlabel(r'$n$')
  plt.ylim(-0.01,0.04)
  plt.xlim(0,0.45)
  #plt.legend(loc=0)

  if nv > 0 and nl > 0:
    # plot a line!
    phi_v = SW.phi(T, nv, npart)
    phi_l = SW.phi(T, nl, npart)
    free_energy_line = (n - nv)*phi_l/(nl - nv) + (n - nl)*phi_v/(nv - nl)
    plt.plot(n*eta_conv, free_energy_line, 'b-')

  if phi_min != 0 or phi_max != 0:
    plt.ylim(phi_min, phi_max)

  plt.savefig(filename)


### Low temp ###

#savefig('figs/SW-phi-lowT.pdf')
plotphi(0.8, 0.04, 'figs/SW-phi-lowT.pdf')

### High Temp ###

#savefig('figs/SW-phi-highT.pdf')
plotphi(1.32, 0.0360850913603, 'figs/SW-phi-highT.pdf')

# Finding mu!

#savefig('figs/fitting-step0-noslope.pdf')
plotphi(0.8005709375, 0.0415310056528, 'figs/fitting-step0-noslope.pdf', 0, 0,
        -0.002, 0.029)

#savefig('figs/fitting-step0.pdf')
plotphi(0.8005709375, 0.0415310056528, 'figs/fitting-step0.pdf', 0.000846555248303, 0.088830194773,
        -0.002, 0.029)

#savefig('figs/fitting-step1-unzoomed.pdf')
plotphi(0.8005709375, 0.0400810433452, 'figs/fitting-step1-unzoomed.pdf', 0.000922047538644, 0.0894990434726,
        -0.002, 0.029)

#savefig('figs/fitting-step1.pdf')
plotphi(0.8005709375, 0.0400810433452, 'figs/fitting-step1.pdf', 0.000922047538644, 0.0894990434726,
        -0.0008, -0.0005)

#savefig('figs/fitting-step2.pdf')
plotphi(0.8005709375, 0.0400859280435, 'figs/fitting-step2.pdf', 0.000921779389983, 0.0894968306258,
        -0.0008, -0.0005)


# T = 0.8005709375
# npart 0.0415310056528
# nv 0.000846555248303
# nl 0.088830194773

# T = 0.8005709375
# npart 0.0400810433452
# nv 0.000922047538644
# nl 0.0894990434726

# T = 0.8005709375
# npart 0.0400859280435
# nv 0.000921779389983
# nl 0.0894968306258
