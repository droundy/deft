#!/usr/bin/python

from __future__ import division

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

from pylab import *
from scipy.special import erf

r0 = 2
R = r0/2
dr = r0/200
rsmall = arange(0, r0/2 + dr/2, dr)
rlarge = arange(r0/2, r0 + dr/2, dr)
r = arange(0,r0 + dr/2, dr)

betaV0 = 100

gamma = betaV0/2.7
alpha = betaV0/r0
A = 4*sqrt(betaV0/r0**2*exp(-betaV0))

#plot(rsmall, A**2*rsmall*exp(alpha*rsmall), 'r-', label="exponential w2 convolved with itself")
#plot(rlarge, A**2*(r0-rlarge)*exp(alpha*rlarge), 'r-')

#plot(r, betaV0/r0*exp(-betaV0*(1-r/r0)), 'b-', label="exponential f'")
gaussianfprime = -2*betaV0*(r-r0)/r0*exp(-betaV0*(1-r/r0)**2)
#gaussianfprime = gaussianfprime/gaussianfprime.max()
plot(r, gaussianfprime, label="gaussian f'")

#plot(rsmall, R*exp(-gamma*(2-2*R/rsmall-(R/rsmall)**2/2))*
#                (2*rsmall*R*sqrt(gamma)*exp(-gamma*rsmall**2/(2*R**2)) +
#                 sqrt(2*pi)*(gamma*rsmall**2-R**2)*erf(sqrt(gamma)/(2*R)*rsmall)))

B = sqrt(1-exp(-betaV0))*2*gamma/(sqrt(pi*gamma)-1)/R**2
alphag = B**2*exp(-gamma*2*(1-rlarge/(2*R))**2)
root2gam = sqrt(2*gamma)
w2w2 = 2*B**2*R/root2gam*(sqrt(pi)/4*(rlarge**2-R**2/gamma)*exp(-2*gamma*(1-rlarge/(2*R))**2)*erf((1-rlarge/(2*R))*root2gam)
                        +R/root2gam*(R-rlarge/2)*exp(-4*gamma*(1-rlarge/(2*R))**2))
#w2w2 = w2w2/w2w2.max()
plot(rlarge, w2w2, label="gaussian W2*W2")
xlim(0,2)
ylim(ymin=0)

# w2w2erfterm = B**2*R/root2gam*(sqrt(pi)/4*(rlarge**2-R**2/gamma)*exp(-2*gamma*(1-rlarge/(2*R))**2)*erf((1-rlarge/(2*R))*root2gam))
# plot(rlarge, w2w2erfterm, label="gaussian W2*W2 erf term")

xlabel('r/R')
legend(loc='best')

savefig('figs/w2convolves.pdf')
