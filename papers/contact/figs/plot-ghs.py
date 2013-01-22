#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
#matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 2:
    print("Usage:  " + sys.argv[0] + " out-filename.pdf")
    exit(1)

pdffilename = sys.argv[1]

pylab.figure(1)
pylab.title('$g_{HS}(r)$')
pylab.axvline(x=1, color='k', linestyle=':')
pylab.axhline(y=1, color='k', linestyle=':')

pylab.figure(2)
pylab.axvline(x=1, color='k', linestyle=':')
pylab.axhline(y=0, color='k', linestyle='-')

def read_ghs(ff):
    mcdatafilename = "figs/mc-inner-sphere-2-0.%d.dat" % (10*ff)
    print 'Using', mcdatafilename
    mcdata = numpy.loadtxt(mcdatafilename)
    nA_mc = mcdata[:,11]
    n0_mc = mcdata[:,10]
    r_mc = mcdata[:,0]
    n_mc = mcdata[:,1]
    ghs = n_mc[r_mc>2]*4*numpy.pi/3/ff
    return r_mc[r_mc>2], ghs

for (c,ff) in [('r',.1008),('b',.2),('g',.3),('c',.4)]:
    r_mc, ghs = read_ghs(ff)

    ghssigma = ghs.max()
    #pylab.plot(r_mc[r_mc>2]/2, (ghs-1)/(ghssigma-1), c+"-",label='filling fraction %.1f'%ff)
    pylab.figure(1)
    pylab.plot(r_mc, ghs, c+"-",label='filling fraction %.1f'%ff)

    ghsnorm = (ghs-1)/(ghssigma-1)
    _, ghs1 = read_ghs(0.1008)
    ghs1sigma = ghs1.max()
    ghs1norm = (ghs1-1)/(ghs1sigma-1)
    #pylab.plot(r_mc, ghs1norm, "r.-")
    #pylab.plot(r_mc, (ghsnorm - ghs1norm)/ghssigma**2, c+"-",label='filling fraction %.1f'%ff)
    # Now a simple model...
    r0 = 0.7
    
    # Below, g1, g2 and g3 are the first three terms in the
    # Percus-Yevik power series for the radial distribution function
    # g(r), taken from the paper "Radial distribution function for
    # hard spheres" by Yuste and Santos from 1991.  It is a power
    # series in *filling fraction*, and the terms are
    # piecewise-defined approximations.
    x = (r_mc)/2 # "diameter units"
    g1 = 8 - 6*x + x**3/2
    g1[x>2] = 0
    #g1 /= g1.max()
    theta1 = 0*x + 1
    theta1[x>=2] = 0
    theta1[x<1] = 0
    theta2 = 0*x + 1
    theta2[x>=3] = 0
    theta2[x<2] = 0
    theta3 = 0*x + 1
    theta3[x>=4] = 0
    theta3[x<3] = 0
    g2 = theta1*(x**6/35 - 9/5*x**4 + 6*x**3 + 9*x**2 - 258*x/5 + 47 - 162/35/x) + \
        theta2*(-x**6/35 + 9/5*x**4 - 6*x**3 - 9*x**2 + 324*x/5 - 81 + 486/35/x)
    r = x
    g3 = theta1*(r**9/2100 - 3/35*r**7 + 16/35*r**6 +12/5*r**5 - 573/25*r**4 + 77/2*r**3 + 2658/35*r**2 - 5367/20*r + 22843/105 - 13299/350/r) +\
         theta2*(-r**9/1050 + 6/35*r**7 - 32/35*r**6 - 24/5*r**5 + 249/5*r**4 - 205/2*r**3 - 1146/7*r**2 + 3423/4*r - 20323/21 + 89049/350/r) +\
         theta3*(1/2100*r**9-3/35*r**7+16/35*r**6+12/5*r**5-672/25*r**4+64*r**3+3072/35*r**2-3072/5*r+90112/105-49152/175/r)
    #pylab.plot(r_mc, 1 + g1*ff + 0*g2*ff**2 + 0*g3*ff**3, c+'--')
    pylab.plot(r_mc, 1 + g1*ff, c+'--')
    
    sigma_correction = (ghs - (1 + g1*ff)).max()
    pylab.figure(2)
    pylab.plot(r_mc, (ghs - (1 + g1*ff)), c+"-",label='filling fraction %.1f'%ff)
    pylab.figure(3)
    pylab.plot(r_mc, (ghs - (1 + g1*ff))/sigma_correction, c+"-",label='filling fraction %.1f'%ff)
    fac = ghssigma-1
    period = 6.0/2 # converting from units of R to units of diameter
    slope = 1.5 # dimensions of 1/distance
    a1 = numpy.exp(-(x-1)*3)*(1 - slope*(period/2/numpy.pi)*fac*numpy.sin(2*numpy.pi*(x-1)/period))
    
    B = 1 - 1/fac
    a1 = (1-B)*numpy.exp(-(x-1)*3) + B*numpy.exp(-(x-1)*6)
    #pylab.plot(r_mc, 1 + a1/a1.max()*(ghssigma-1), c+"-.")
    fac = (ghs - (1 + g1*ff)).max()
    slope = fac/0.5
    #pylab.plot(r_mc, fac*numpy.exp(-(x-1)*2/0.25), c+"-.")
    correction = fac*(1 - (x-1)*2/0.35)*numpy.exp(-((x-1)*2/0.5))
    correction = fac*(1 - (x-1)*2/0.35)*numpy.exp(-((x-1)*2/0.8)**2)
    minloc = 3.0
    minval = -0.5
    correction = fac*( (1-minval)*(2*x - minloc)**2/(minloc-2)**2 + minval )
    maxr = 4.0
    correction[r_mc>maxr] = 0
    # The following is motivated by the paper of Trokhymchuk et al,
    # and uses their notation.
    sigma = 2 # the diameter
    d = (2*ff)**(1/3) # ???
    mu = -d/sigma
    beta0 = - numpy.sqrt(3)*d/2/sigma
    alpha0 = d/2/sigma
    gamma = -sigma/beta0*((alpha0*sigma-mu*sigma)*(alpha0**2+beta0**2) +
                          alpha0**2 + beta0**2 - mu*alpha0)
    omega = (4.720 - 0.682)/sigma
    kappa = (4.67 + 3.536)/sigma
    alpha = 44.554
    beta = -5.022
    rstar = sigma*2.0016
    gm = 1.0286
    delta = -omega*rstar - (kappa*rstar + 1)/(omega*rstar)
    B = rstar*(gm - ghssigma/rstar*numpy.exp(mu*(rstar-sigma)))/  \
        (numpy.cos(beta*(rstar-sigma)+gamma)*numpy.exp(alpha*(rstar-sigma)) - \
         numpy.cos(gamma)*numpy.exp(mu*(rstar-sigma)))
    A = sigma*ghssigma - B*numpy.cos(gamma)
    C = rstar*(gm-1)*numpy.exp(kappa*rstar)/numpy.cos(omega*rstar+delta)
    correction = A/r_mc*numpy.exp(mu*(r_mc-sigma)) + \
                 B/r_mc*numpy.cos(beta*(r_mc-sigma)+gamma)*numpy.exp(alpha*(r_mc-sigma))
    correction2 = 1 + C/r_mc*numpy.cos(omega*r_mc+delta)*numpy.exp(-kappa*r_mc)
    correction[r_mc>rstar] = correction2[r_mc>rstar]
    correction[r_mc<sigma] = 0
    # Now here is my very hokey function
    #pylab.figure(1)
    slope = -3.4
    mu = -1.7
    beta = 2*numpy.pi/(2.2*sigma)
    alpha = -.8
    frac = (slope - mu)/beta
    gamma = 0
    correction = numpy.exp(mu*(r_mc-sigma))*sigma/r_mc + \
        frac*(beta*(r_mc-sigma))*numpy.exp(-(alpha*(r_mc-sigma))**2)*sigma/r_mc
    correction[r_mc<sigma] = 0
    correction /= correction.max()
    # Another hokey attempt
    A = 4.0
    a = 0.5
    B = -4.0
    b = 1.0
    C = 1 - A - B
    cc = 1.5
    correction = A*numpy.exp(-(r_mc-sigma)/a) +  B*numpy.exp(-(r_mc-sigma)/b) +  C*numpy.exp(-(r_mc-sigma)/cc)
    correction[r_mc<sigma] = 0
    correction /= correction.max()
    # Another hokey attempt
    A = 1.0
    a = 0.8
    B = sigma*(slope + A/a)
    b = 0.6
    correction = A*numpy.exp(-(r_mc-sigma)/a) +  B*numpy.exp(-(r_mc-sigma)/b)*(r_mc-sigma)/sigma
    correction[r_mc<sigma] = 0
    correction /= correction.max()
    # Another hokey attempt
    a = 0.6
    b = 4.0
    correction = numpy.exp(-(r_mc-sigma)/a)*(1 + slope*(r_mc-sigma) + b*((r_mc-sigma)/sigma)**2)
    correction[r_mc<sigma] = 0
    correction /= correction.max()
    pylab.figure(2)
    pylab.plot(r_mc, correction*(ghssigma-1-g1.max()*ff), c+"-.")
    pylab.figure(3)
    pylab.plot(r_mc, correction, c+"-.")

print r_mc[0]

pylab.figure(1)
pylab.xlim(2,6.5)
pylab.xlabel(r"$r/R$")
pylab.ylabel("$g(r)$")
pylab.legend(loc='best').get_frame().set_alpha(0.5)
pylab.savefig(pdffilename)

pylab.figure(2)
pylab.xlim(2,6.5)
pylab.xlabel(r"$r/R$")
pylab.ylabel("$g(r)$")
pylab.legend(loc='best').get_frame().set_alpha(0.5)

pylab.show()
