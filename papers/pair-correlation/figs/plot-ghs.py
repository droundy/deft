#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys

if len(sys.argv) < 2 or sys.argv[1]!="show" :
  matplotlib.use('Agg')
import pylab, numpy, random
from pylab import *
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
x[0] = 7.9

x[1] = -.13
x[2] = 7.13/2
x[3] = 1.71

x[4] = .18
x[5] = 5.52
x[6] = 3.39

# x[7] = .03
# x[8] = 4.56
# x[9] = 3.32

colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5]
able_to_read_file = True


def read_ghs(base, ff):
    mcdatafilename = "%s-0.%02d.dat" % (base, 100*ff)
    try:
        mcdata = numpy.loadtxt(mcdatafilename)
    except IOError:
        global able_to_read_file
        able_to_read_file = False
        return 0,0
    print 'Using', mcdatafilename, 'for filling fraction', ff
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
    if able_to_read_file == False:
        break
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

if able_to_read_file == False:
    matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
    matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
    pylab.savefig("figs/ghs-g.pdf")
    pylab.savefig("figs/ghs-g-ghs.pdf")
    exit(0)

def evalg(x, gsigma, r):
  hsigma_rolloff = 5.0
  hsigma = gsigma - 1 # (1 - 0.5*eta)/(1-eta)**3 - 1
  h0 = hsigma # was x[0]*gsig
  f0 = numpy.exp(-x[0]*r)

  # the slope dghs/dr should be = -hsigma - 0.7*hsigma**2 according my fit (by eye)
  # but our "r" is r/2, so our slope is twice that.
  # -2*hsigma - 1.4*hsigma**2 = -x[0]*hsigma -amplitude[4]*x[5]
  # thus amplitude[4] = (hsigma*(2-x[0]) + 1.4*hsigma**2)/x[5]

  h1 = x[1]*hsigma
  #h1 = x[1]*hsigma_rolloff*(1-exp(-hsigma/hsigma_rolloff))
  f1 = numpy.sin(x[2]*r)**2 * numpy.exp(-x[3]*r)
  h2 = -x[4]*(hsigma + hsigma**2)
  h2 = (hsigma*(-2+x[0]) -2*hsigma**2)/x[5]
  #h2 = -x[4]*hsigma_rolloff**2*(1-exp(-hsigma**2/hsigma_rolloff**2))
  f2 = numpy.sin(x[5]*r) * numpy.exp(-x[6]*r)
  #h3 = -x[7]*hsigma**(3)
  #f3 = numpy.sin(x[8]*r) * numpy.exp(-x[9]*r)
  return 1 + h0*f0 + h1*f1 + h2*f2#+ h3*f3

def dist(x):
    # function with x[i] as constants to be determined
    g = numpy.zeros_like(gsigconcatenated)
    for i in range(len(g)):
      g[i] = evalg(x, gsigconcatenated[i], rconcatenated[i])
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
#vals = x
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

pylab.figure()
for eta in [.5, .6, .7, .8]:
  gsigma = (1-eta/2)/(1-eta)**3
  pylab.plot(r_mc, evalg(vals, gsigma, r), label='eta %g  gsig %g'%(eta, gsigma))
axhline(y=0)
pylab.xlim(2,6.5)
pylab.legend(loc='best')
pylab.show()



