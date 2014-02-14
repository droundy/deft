#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
from pylab import *

from scipy.special import erf

diameter = 2
R = 1
V0 = 1

dr = R/100
r = arange(0,1.5*diameter, dr)


V = V0*(1-r/diameter)**2
V[r>diameter] = 0

figure()
plot(r, V)
title('$V(r)$')

kT = 0.2*V0
beta = 1/kT
f = exp(-beta*V) - 1


sigma = diameter*(1 - sqrt(kT/V0*log(2)))
#Vprime_sigma = 2*V0*(1-sigma/diameter)*(-1/diameter)
delta_r = diameter*sqrt(kT/V0/log(2)/pi)
C = 0.5

figure()
plot(r,f)
ferf = C*(erf((r-sigma)/delta_r)-1)
plot(r, ferf)
axvline(sigma)
title('$f(r)$')

Vprime = 2*V0*(1-r/diameter)*(-1/diameter)
fprime = -beta*Vprime*exp(-beta*V)
fprime[r>diameter] = 0

figure()
plot(r,fprime, 'b-', label="true $f'(r)$ for harmonic potential")
plot(r, C*exp(-(r - sigma)**2/delta_r**2)/delta_r*2/sqrt(pi), 'g-', label='erf approximation')

axvline(sigma, color='k', linestyle=':')
xlim(0, 1.1*diameter)
legend(loc='best').draw_frame(False)
title("$f'(r)$")
savefig('figs/erf.pdf')

figure()
plot(r, -kT*log(0.5*(erf((r-sigma)/delta_r)+1)), label='erf approximation')
plot(r[r<2*R], V0*(1 - r[r<2*R]/(2*R))**2, label='harmonic potential')
axvline(sigma, color='k', linestyle=':')
axhline(kT*log(2), color='k', linestyle=':')
ylim(0, 4*kT)
savefig('figs/erf-potential.pdf')

show()
