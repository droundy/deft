#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, random
#import Scientific.Functions.LeastSquares as ls

from scipy.optimize import leastsq

pylab.figure(1)
pylab.title('$g_{HS}(r)$') #radial distribution for hard spheres
pylab.axvline(x=1, color='k', linestyle=':')
pylab.axhline(y=1, color='k', linestyle=':')

#pylab.figure(2)
#pylab.axvline(x=1, color='k', linestyle=':')
#pylab.axhline(y=0, color='k', linestyle='-')

toFit = 4
numParams = 7+3
x = [1]*numParams
# x[0] = 1.0 # 1.15
x[0] = 6.134

x[1] = .413
x[2] = 6.449
x[3] = 2.058

x[4] = .331
x[5] = 4.560
x[6] = 3.319

x[7] = .0331
x[8] = 4.560
x[9] = 3.319

colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5]



def read_ghs(base, ff):
    mcdatafilename = "%s-0.%02d.dat" % (base, 100*ff)
    print 'Using', mcdatafilename, 'for filling fraction', ff
    mcdata = numpy.loadtxt(mcdatafilename)
    r_mc = mcdata[:,0]
    n_mc = mcdata[:,1]
    ghs = n_mc/ff
    return r_mc, ghs


# READ DATA
ghs = [0]*len(ff)
#ghslores = [0]*len(ff)
gsig = [0]*len(ff)
i = 0
while (i < len(ff)):
    r_mc, ghs[i] = read_ghs("figs/gr", ff[i])
    #r_mclores, ghslores[i] = read_ghs("grlores", ff[i])
    pylab.figure(1)
    pylab.plot(r_mc, ghs[i], colors[i]+"-",label='ghs at filling fraction %.2f'%ff[i])
    # The following is the Monte Carlo approximation of the
    # distribution function at contact.  This gives us an answer with
    # no systematic error (well, very little, and what we have is due
    # to the finite bin size), and if we wait long enough the
    # statistical error should become negligible.
    #gsig[i] = ghs[i].max()
    # The following is the Carnahan Starling formula for the
    # distribution function at contact.  This is approximate, but it
    # is a good approximation.
    gsig[i] = (1 - 0.5*ff[i])/(1-ff[i])**3
    r = (r_mc)/2 - 1 # "diameter units", shifted to x = 0 at d = 1
    i += 1

def dist(x):
    # function with x[i] as constants to be determined
    g = numpy.zeros_like(gsigconcatenated)
    for i in range(len(g)):
        # note: using rconcatenated = r-2, so it returns exactly gsigma at contact
        hsigma = gsigconcatenated[i] - 1
        h0 = hsigma # was x[0]*gsig
        f0 = numpy.exp(-x[0]*rconcatenated[i])
        h1 = x[1]*hsigma
        f1 = numpy.sin(x[2]*rconcatenated[i]) * numpy.exp(-x[3]*rconcatenated[i])
        h2 = -x[4]*hsigma**(2)
        f2 = numpy.sin(x[5]*rconcatenated[i]) * numpy.exp(-x[6]*rconcatenated[i])
        h3 = -x[7]*hsigma**(3)
        f3 = numpy.sin(x[8]*rconcatenated[i]) * numpy.exp(-x[9]*rconcatenated[i])
        g[i] = 1 + h0*f0 + 0*h1*f1 + 0*h2*f2 + 0*h3*f3
    return g

def dist2(x):
    return dist(x) - ghsconcatenated

ghsconcatenated = ghs[0]
for i in range(1,len(ff)):
    ghsconcatenated = numpy.concatenate((ghsconcatenated, ghs[i]))

gsigconcatenated = [0]*len(r)*len(gsig)
j = 0
while (j < len(gsig)):   # makes ind an array of pairs, ind = [(sig1, r1), (sig1, r2), ..., (sig2, r2), ...]
    i = 0
    while (i < len(r)):
        gsigconcatenated[i + j*len(r)] = gsig[j]
        i += 1
    j += 1

rconcatenated = [0]*len(r)*len(gsig)
j = 0
while (j < len(gsig)):   # makes ind an array of pairs, ind = [(sig1, r1), (sig1, r2), ..., (sig2, r2), ...]
    i = 0
    while (i < len(r)):
        rconcatenated[i + j*len(r)] = r[i]
        i += 1
    j += 1



vals = [1]*numParams

print "beginning least squares fit..."
vals, cov = leastsq(dist2, x)
print "original fit complete, cov: " + str(cov)

i = 0
while (i < len(vals)):
    print 'vals[' + str(i) + ']: ' + str(vals[i])
    i += 1

for i in numpy.arange(len(vals)):
    vals[i] = round(vals[i],2)
i=0
while (i < len(vals)):
    print 'x[' + str(i) + ']: ' + str(vals[i])
    i += 1

g = dist(vals)
gdifference = dist2(vals)
# for i in g:
#     print dist(vals, ind)[i]

for i in range(len(ff)):
    pylab.figure(1)
    pylab.plot(r_mc, g[i*len(r):(i+1)*len(r)], colors[i]+'--',label='g at filling fraction %.2f'%ff[i])

    pylab.figure(2)
    pylab.plot(r_mc, gdifference[i*len(r):(i+1)*len(r)], colors[i]+'--')
    pylab.plot(r_mc, g[i*len(r):(i+1)*len(r)] - ghs[i], colors[i]+'-')
    #pylab.plot(r_mc, numpy.abs(numpy.asarray(ghsconcatenated[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'-')



pylab.figure(1)
pylab.xlim(2,6.5)
#pylab.ylim(0.5, 3.5)
pylab.xlabel(r"$r/R$")
pylab.ylabel("$g(r)$")
pylab.legend(loc='best').get_frame().set_alpha(0.5)

pylab.savefig("figs/ghs-g.pdf")


pylab.figure(2)
pylab.xlim(2,6.5)
#pylab.ylim(0, 2.0)
pylab.xlabel(r"$r/R$")
pylab.ylabel("|ghs - g|")
#pylab.legend(loc='best').get_frame().set_alpha(0.5)
pylab.savefig("figs/ghs-g-ghs.pdf")

pylab.show()



