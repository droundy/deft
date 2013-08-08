#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys, os.path

if len(sys.argv) < 2 or sys.argv[1]!="show" :
  matplotlib.use('Agg')

from pylab import *
from scipy.optimize import leastsq

figure(1)
title('$g_{HS}(r)$') #radial distribution for hard spheres
axvline(x=1, color='k', linestyle=':')
axhline(y=1, color='k', linestyle=':')

x = zeros(7)

x[1] = -.13
x[2] = 7.13/2
x[3] = 1.71

x[1] = 0.435
x[2] = 3.552
x[3] = 0.870

x[4] = 0.329
x[5] = 2.540
x[6] = 1.940

colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = array([.05, .1, .15, .2, .25, .3, .35, .4, .45, .5])
able_to_read_file = True


def read_ghs(base, ff):
    mcdatafilename = "%s-%4.2f.dat" % (base, ff)
    if (os.path.isfile(mcdatafilename) == False):
      print "File does not exist: ", mcdatafilename
      global able_to_read_file
      able_to_read_file = False
      return 0, 0

    mcdata = loadtxt(mcdatafilename)
    print 'Using', mcdatafilename, 'for filling fraction', ff
    r_mc = mcdata[:,0]
    n_mc = mcdata[:,1]
    ghs = n_mc/ff
    return r_mc, ghs


# READ DATA
ghs = [0]*len(ff)
#ghslores = [0]*len(ff)
eta = [0]*len(ff)

for i in range(len(ff)):
    r_mc, ghs[i] = read_ghs("figs/gr", ff[i])
    if able_to_read_file == False:
        break
    #r_mclores, ghslores[i] = read_ghs("grlores", ff[i])
    figure(1)
    plot(r_mc, ghs[i], colors[i]+"-",label='ghs at filling fraction %.2f'%ff[i])
    # The following is the Monte Carlo approximation of the
    # distribution function at contact.  This gives us an answer with
    # no systematic error (well, very little, and what we have is due
    # to the finite bin size), and if we wait long enough the
    # statistical error should become negligible.
    #gsig[i] = ghs[i].max()
    # The following is the Carnahan Starling formula for the
    # distribution function at contact.  This is approximate, but it
    # is a good approximation.
    eta[i] = ff[i]
    r = r_mc - 2 # shift to x = 0 at contact

if able_to_read_file == False:
    plot(arange(0,10,1), [0]*10, 'k')
    suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
    savefig("figs/ghs-g.pdf")
    savefig("figs/ghs-g-ghs.pdf")
    exit(0)

def evalg(x, eta, r):
  hsigma_rolloff = 5.0
  hsigma = (1 - 0.5*eta)/(1-eta)**3 - 1
  h0 = hsigma # was x[0]*gsig

  f0 = exp(-x[0]*r)

  # the slope dghs/dr should be = -hsigma - 0.7*hsigma**2 according my fit (by eye)
  # but our "r" is r/2, so our slope is twice that.
  # -2*hsigma - 1.4*hsigma**2 = -x[0]*hsigma -amplitude[4]*x[5]
  # thus amplitude[4] = (hsigma*(2-x[0]) + 1.4*hsigma**2)/x[5]

  h1 = x[1]*hsigma
  #h1 = x[1]*hsigma_rolloff*(1-exp(-hsigma/hsigma_rolloff))
  f1 = sin(x[2]*r)**2 * exp(-x[3]*r)
  h2 = -x[4]*(hsigma + hsigma**2)
  h2 = (hsigma*(-2+x[0]) -2*hsigma**2)/x[5]
  #h2 = -x[4]*hsigma_rolloff**2*(1-exp(-hsigma**2/hsigma_rolloff**2))
  f2 = sin(x[5]*r) * exp(-x[6]*r)

  #h3 = -x[7]*hsigma**(3)
  #f3 = sin(x[8]*r) * exp(-x[9]*r)
  return 1 + h0*f0 + h1*f1 + h2*f2 #+ h3*f3

def dist(x):
    # function with x[i] as constants to be determined
    g = zeros_like(etaconcatenated)
    for i in range(len(g)):
      g[i] = evalg(x, etaconcatenated[i], rconcatenated[i])
    return g

def dist2(x):
    return dist(x) - ghsconcatenated


ghsconcatenated = ghs[0]
for i in range(1,len(ff)):
    ghsconcatenated = concatenate((ghsconcatenated, ghs[i]))

etaconcatenated = [0]*len(r)*len(eta)
j = 0
while (j < len(eta)):
    i = 0
    while (i < len(r)):
        etaconcatenated[i + j*len(r)] = eta[j]
        i += 1
    j += 1

rconcatenated = [0]*len(r)*len(eta)
j = 0
while (j < len(eta)):
    i = 0
    while (i < len(r)):
        rconcatenated[i + j*len(r)] = r[i]
        i += 1
    j += 1

vals = zeros_like(x)

print "beginning least squares fit..."
vals, mesg = leastsq(dist2, x)
chi2 = sum(dist2(vals)**2)
print "original fit complete, chi^2: %.3f" % chi2


for i in range(len(x)):
  print "vals[%i]: %.3f\t x[%i]: %g" %(i, vals[i], i, x[i])

g = dist(vals)
gdifference = dist2(vals)

# for i in g:
#     print dist(vals, ind)[i]

for i in range(len(ff)):
    figure(1)
    plot(r_mc, g[i*len(r):(i+1)*len(r)], colors[i]+'--',label='g at filling fraction %.2f'%ff[i])
    figure(2)
    plot(r_mc, gdifference[i*len(r):(i+1)*len(r)], colors[i]+'--')
    plot(r_mc, g[i*len(r):(i+1)*len(r)] - ghs[i], colors[i]+'-')
    #plot(r_mc, numpy.abs(numpy.asarray(ghsconcatenated[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'-')



figure(1)
xlim(2,6.5)
#ylim(0.5, 3.5)
xlabel(r"$r/R$")
ylabel("$g(r)$")
legend(loc='best').get_frame().set_alpha(0.5)

savefig("figs/ghs-g.pdf")


figure(2)
xlim(2,6.5)
#ylim(0, 2.0)
xlabel(r"$r/R$")
ylabel("|ghs - g|")
#legend(loc='best').get_frame().set_alpha(0.5)
savefig("figs/ghs-g-ghs.pdf")

figure(3)
for eta in [.5, .6, .7, .8]:
  gsigma = (1-eta/2)/(1-eta)**3
  plot(r_mc, evalg(vals, eta, r), label='eta %g  gsig %g'%(eta, gsigma))
axhline(y=0)
xlim(2,6.5)
legend(loc='best')
show()



