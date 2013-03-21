#!/usr/bin/python

from __future__ import division

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

from pylab import *
from scipy.special import erf
from itertools import *
from scipy.special import erf

V0 = 1
Temp = 0.01
betaV0 = V0/Temp
r0=1

dr = 0.0001
r = arange(0, 1.2 + dr/2, dr)

def findgamma(betaV0):
    return 2*((sqrt(pi*betaV0)+sqrt(pi*betaV0-16*sqrt(betaV0)))/8)**2

gamma = 2*((sqrt(pi*betaV0)+sqrt(pi*betaV0-16*sqrt(betaV0)))/8)**2
sg = sqrt(gamma)

w3 = -1/(sqrt(pi*gamma) - 1)*(1 - exp(-gamma*(1-r)**2) - sqrt(pi*gamma)*erf(sg*(1-r)))
def w3(betaV0):
    gamma = findgamma(betaV0)
    sg = sqrt(gamma)
    return -1/(sqrt(pi*gamma) - 1)*(1 - exp(-gamma*(1-r)**2) - sqrt(pi*gamma)*erf(sg*(1-r)))

figure(figsize=(4.5,4))
plot(r, w3(betaV0), 'b-', label='$kT=%g V_{max}$' % (1/betaV0))
plot(r, w3(betaV0*10), 'r-', label='$kT=%g V_{max}$' % (0.1/betaV0))
ylim(0,1.1)
legend(loc='best')
ylabel('$w_3(r)$')
savefig('figs/w_3.pdf')
clf()

w2 = (2*gamma*r)/((sqrt(pi*gamma)-1)*r0)*exp(-gamma*(1-r)**2)

def w2(betaV0):
    gamma = findgamma(betaV0)
    sq = sqrt(gamma)
    func = (2*gamma*r)/((sqrt(pi*gamma)-1)*r0)*exp(-gamma*(1-r)**2)
    func[r>1] = 0
    return func

figure(figsize=(4.5,4))
plot(r, w2(betaV0), 'b-', label='$kT=%g V_{max}$' % (1/betaV0))
plot(r, w2(betaV0*10), 'r-', label='$kT=%g V_{max}$' % (0.1/betaV0))
ylim(ymin=0)
legend(loc='best')
ylabel('$w_2(r)$')
savefig('figs/w_2.pdf')
clf()

w1 = (gamma/r)/(2*(sqrt(pi*gamma)-1)*r0**2)*exp(-gamma*(1-r)**2)

def w1(betaV0):
    gamma = findgamma(betaV0)
    sq = sqrt(gamma)
    func = (gamma/r)/(2*(sqrt(pi*gamma)-1)*r0**2)*exp(-gamma*(1-r)**2)
    func[r>1] = 0
    return func

figure(figsize=(4.5,4))
plot(r, w1(betaV0), 'b-', label='$kT=%g V_{max}$' % (1/betaV0))
plot(r, w1(betaV0*10), 'r-', label='$kT=%g V_{max}$' % (0.1/betaV0))
ylim(ymin=0)
legend(loc='best')
ylabel('$w_0(r)$')
savefig('figs/w_1.pdf')
clf()

stepfunc = 0*r
stepfunc[r<=1] = 1

plot(r,stepfunc)
ylim(0,1.1)
ylabel('$\Theta(R-r)$')
xlabel('$r/R$')
savefig('figs/step.pdf')
clf()

stepfunc *= 0
stepfunc[r == 1] = 1e300
plot(r,stepfunc)
ylabel('$\delta(R-r)$')
xlabel('$r/R$')
savefig('figs/delta.pdf')
clf()
