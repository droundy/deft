#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys, os.path

if len(sys.argv) < 2 or sys.argv[1]!="show" :
  matplotlib.use('Agg')

from pylab import *
from scipy.optimize import leastsq

sigma = 2

figure(1)
title('$g_{HS}(r)$') #radial distribution for hard spheres
axvline(x=sigma, color='k', linestyle=':')
axhline(y=1, color='k', linestyle=':')

x = zeros(6)

x[0] = 1.5
x[1] = 1.5
x[2] = 1.6
x[4] = 1e-5
x[5] = 0.8

x[0] = 1.706
x[1] = 1.358
x[2] = 1.567
x[4] = 0.675
x[5] = 0.463

x[0] = 1.854
x[1] = 1.436
x[2] = 1.445
x[4] = 0.424
x[5] = 0.330

x[0] = 2.033
x[1] = 1.584
x[2] = 1.503
x[4] = 0.329
x[5] = 0.269

x[0] = 2.143
x[1] = 1.588
x[2] = 1.422
x[4] = 0.268
x[5] = 0.227

colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = array([0.05, .1, 0.15, .2, .25, .3, .35,  .4, .45])
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
    plot(r_mc, ghs[i], colors[i]+":",label='$\eta = %.2f$'%ff[i])
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
    r = r_mc # shift to x = 0 at contact

if able_to_read_file == False:
  plot(arange(0,10,1), [0]*10, 'k')
  suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
  savefig("figs/ghs-g3.pdf")
  savefig("figs/ghs-g-ghs3.pdf")
  exit(1)

toprint = True
def evalg2(xnew, eta, r):
  global toprint
  z = r - sigma
  hsigma = (1 - 0.5*eta)/(1-eta)**3 - 1
  density = 3/4/pi*eta
  rhs = (1-eta)**4/(1+4*eta+4*eta**2-4*eta**3+eta**4)/3

  x0 = xnew[0]
  x1 = xnew[1]
  x2 = xnew[2]
  x4 = xnew[4]
  x5 = xnew[5]

  int_h0 = 4*pi*hsigma*(2 + sigma*x0*(2 + sigma*x0))/x0**3

  int_h1_over_A = 4*pi*x1* \
      (sigma**2*x1**4 + \
         2*x1**2*(-1 + sigma*x2*(2 + sigma*x2)) + \
         x2**2*(6 + sigma*x2*(4 + sigma*x2)) \
         )/(x1**2 + x2**2)**3

  int_h2_over_B = 4*pi*x4* \
      (sigma**2*x4**4 + 2*x4**2* \
         (-1 + sigma*x5*(2 + sigma*x5)) + x5**2*(6 + sigma*x5*(4 + sigma*x5))) \
         /(x4**2 + x5**2)**3

  A = ((rhs-1)/density - int_h0 - ((x0-1)*hsigma-hsigma**2)/x4*int_h2_over_B) \
      / (int_h1_over_A - x1/x4*int_h2_over_B)
  B = ((x0-1)*hsigma - hsigma**2)/x4 - A*x1/x4

  if toprint:
    print eta, A, B
    toprint = False

  f0 = hsigma
  j0 = exp(-x[0]*z)

  j1 = A*sin(xnew[1]*z)*exp(-x[2]*z)

  j2 = B*sin(xnew[4]*z)*exp(-x[5]*z)

  #f3 = -xnew[6]*hsigma**3
  #j3 = sin(xnew[7]*r)*exp(-xnew[8]*r)

  return 1 + f0*j0 + j1 + j2# + f3*j3

def dist(x):
  # function with x[i] as constants to be determined
  g = zeros_like(etaconcatenated)
  for i in range(len(g)):
    g[i] = evalg2(x, etaconcatenated[i], rconcatenated[i])
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

chi2initially = sum(dist2(x)**2)
print "chi^2 (initial guess, initially) =", chi2initially

vals = zeros_like(x)

if True:
  print "beginning least squares fit..."
  vals, cov = leastsq(dist2, x)
  print "original fit complete, cov: " + str(cov)
else:
  print "skipping the least squares fit..."
  vals = x

toprint = True
for i in range(len(x)):
  print "vals[%i]: %.3f\t x[%i]: %g" %(i, vals[i], i, x[i])

g = dist(vals)
valsrounded = zeros_like(vals)
for i in range(len(valsrounded)):
  valsrounded[i] = float("%.3f" % vals[i])
grounded = dist(valsrounded)
gdifference = dist2(vals)

print "chi^2 (initial guess) =", sum(dist2(x)**2)
print "chi^2 =", sum(dist2(vals)**2)
toprint = True
chi2rounded = sum(dist2(valsrounded)**2)
print "chi^2 (rounded) =", chi2rounded

# for i in g:
#     print dist(vals, ind)[i]

for i in range(len(ff)):
  figure(1)
  #plot(r_mc, g[i*len(r):(i+1)*len(r)], colors[i]+'--',label='g at filling fraction %.2f'%ff[i])
  plot(r_mc, grounded[i*len(r):(i+1)*len(r)], colors[i]+'-')
  hsigma = (1 - 0.5*ff[i])/(1-ff[i])**3 - 1
  density = 4/3*pi*ff[i]
  rhs = (1-ff[i])**4/(1+4*ff[i]+4*ff[i]**2-4*ff[i]**3+ff[i]**4)/3
  #integral = hsigma*(1/a + x[0]*x[1]/())



  #print density, integral, rhs
  #print "ff: %.2f\t thing: %g" %(ff[i], 1 - rho*integral - rhs)
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

#axvline(sigma + pi/vals[1], color='r')
#axvline(sigma + pi/vals[4], color='b')

#axvline(sigma + 2*pi/vals[1], color='r')
#axvline(sigma + 2*pi/vals[4], color='b')

savefig("figs/ghs-g3.pdf")


figure(2)
xlim(2,6.5)
#ylim(0, 2.0)
xlabel(r"$r/R$")
ylabel("|ghs - g|")
#legend(loc='best').get_frame().set_alpha(0.5)
savefig("figs/ghs-g-ghs3.pdf")

# figure(3)
# for eta in [.5, .6, .7, .8]:
#   gsigma = (1-eta/2)/(1-eta)**3
#   plot(r_mc, evalg2(vals, eta, r), label='eta %g  gsig %g'%(eta, gsigma))
axhline(y=0)
xlim(2,6.5)
legend(loc='best')
show()



