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

x = zeros(3)

x[0] = 1.220
x[1] = 2.125
x[2] = 2.132

sigma = 2.0
colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = array([.05, .1, .15, .2, .25, .3, .35, .4, .45, .5])
ff = array([.1, .2, .3, .4])
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
    plot(r_mc, ghs[i], colors[i]+":")
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
    r = r_mc

if able_to_read_file == False:
    plot(arange(0,10,1), [0]*10, 'k')
    suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
    savefig("figs/ghs-g.pdf")
    savefig("figs/ghs-g-ghs.pdf")
    exit(0)

def evalg(xnew, eta, r):
  z = r - sigma
  hsigma = (1 - 0.5*eta)/(1 - eta)**3 - 1
  density = 3/4/pi*eta
  rhs = (1-eta)**4/(1 + 4*eta + 4*eta**2 - 4*eta**3 + eta**4)

  a0 = xnew[0]
  a1 = 0
  a2 = xnew[1]
  a3 = 0
  a4 = xnew[2]

  int_h0 = 4*pi*hsigma*(2 + sigma*a0*(2 + sigma*a0))/a0**3 - 4/3*pi*sigma**3

  int_h1_over_a1 = 4*pi*(6 + sigma*a2*(4 + sigma*a2))/a2**4

  int_h2_over_a3 = 8*pi*(12 + sigma*a4*(6 + sigma*a4))/a4**5

  A = (rhs-1)/density - int_h0
  B = int_h2_over_a3
  C = int_h1_over_a1

  a1 = hsigma*(a0 - 1 - hsigma) # sets slope at sigma

  a3 = (A - C*a1)/B # sets integral

  f0 = hsigma
  j0 = exp(-a0*z)

  f1 = a1
  j1 = z*exp(-a2*z)

  f2 = a3
  j2 = z**2*exp(-a4*z)

  return 1 + f0*j0 + f1*j1 + f2*j2

def dist(x):
  # function with x[i] as constants to be determined
  R, ETA = meshgrid(r, eta)
  g = zeros_like(ETA)
  g = evalg(x, ETA, R)
  return reshape(g, len(eta)*len(r))

def dist2(x):
  return dist(x) - reshape(ghs, len(eta)*len(r))


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
    plot(r_mc, g[i*len(r):(i+1)*len(r)], colors[i]+'-',label='$\eta=%.1f$'%ff[i])
    figure(2)
    plot(r_mc, gdifference[i*len(r):(i+1)*len(r)], colors[i]+'--')
    plot(r_mc, g[i*len(r):(i+1)*len(r)] - ghs[i], colors[i]+'-')
    #plot(r_mc, numpy.abs(numpy.asarray(ghsconcatenated[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'-')



figure(1)
xlim(2,5)
#ylim(0.5, 3.5)
xlabel(r"$r/R$")
ylabel("$g(r)$")
legend(loc='best')
#legend(loc='best').get_frame().set_alpha(0.5)

savefig("figs/ghs-g.pdf")


figure(2)
xlim(2,6.5)
#ylim(0, 2.0)
xlabel(r"$r/R$")
ylabel("|ghs - g|")
#legend(loc='best').get_frame().set_alpha(0.5)
savefig("figs/ghs-g-ghs.pdf")

# figure(3)
# for eta in [.5, .6, .7, .8]:
#   gsigma = (1-eta/2)/(1-eta)**3
#   plot(r_mc, evalg(vals, eta, r), label='eta %g  gsig %g'%(eta, gsigma))
# axhline(y=0)
# xlim(2,6.5)
show()



