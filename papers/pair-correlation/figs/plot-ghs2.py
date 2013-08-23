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


def evalg(xnew, eta, r):
  return evalg_quadratic(xnew, eta, r)

x = rand(3)*5

# # Linear:
# x[0] = 2.079
# x[1] = -1.234
# x[2] = 1.990
# x[3] = 1.705
# x[4] = 0

# Quadratic:
x[0] = 0.609
x[1] = 0.978
x[2] = 1.344


colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']*2
ff = array([.1, .2, .3, .4])
#ff = array([.05, .1, .15, .2, .25, .3, .35, .4, .45, .5])
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
    r = r_mc

if able_to_read_file == False:
  plot(arange(0,10,1), [0]*10, 'k')
  suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
  savefig("figs/ghs-g2.pdf")
  savefig("figs/ghs-g-ghs2.pdf")
  exit(0)

def evalg_quadratic(xnew, eta, r):
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


def evalg_nosquare(xnew, eta, r):
  global toprint
  z = r - sigma
  hsigma = (1 - 0.5*eta)/(1-eta)**3 - 1
  density = 3/4/pi*eta
  rhs = (1-eta)**4/(1+4*eta+4*eta**2-4*eta**3+eta**4)

  a0 = xnew[0]
  a1 = 0 # this parameter is constrained via integral
  a2 = xnew[1]
  a3 = xnew[2]
  a4 = 0
  a5 = xnew[3]
  a6 = xnew[4]

  int_h0 = 4*pi*hsigma*(2 + sigma*a0*(2 + sigma*a0))/a0**3 - 4/3*pi*sigma**3

  int_h1_over_a1 = 4*pi*a2* \
      (sigma**2*a2**4 + 2*a2**2* \
         (-1 + sigma*a3*(2 + sigma*a3)) + a3**2*(6 + sigma*a3*(4 + sigma*a3))) \
         /(a2**2 + a3**2)**3

  int_h2_over_a4 = 4*pi*a5* \
      (sigma**2*a5**4 + 2*a5**2* \
         (-1 + sigma*a6*(2 + sigma*a6)) + a6**2*(6 + sigma*a6*(4 + sigma*a6))) \
         /(a5**2 + a6**2)**3

  A = (rhs-1)/density - int_h0
  B = int_h2_over_a4
  C = int_h1_over_a1

  a1 = (A/C - B*hsigma/C/a5*(-1-hsigma + a0))/(1 - B*a2/C/a5)

  a4 = hsigma/a5*(-1 - hsigma + a0 - a1*a2/hsigma) # this matches slope at x = sigma


  if toprint:
    print "a1: %.3f, a4: %.3f" %(a1, a4)
    toprint = False

  f0 = hsigma
  j0 = exp(-a0*z)

  f1 = a1
  j1 = sin(a2*z)*exp(-a3*z)

  f2 = a4
  j2 = sin(a5*z)*exp(-a6*z)

  return 1 + f0*j0 + f1*j1 + f2*j2

def evalg_linear(xnew, eta, r):
  global toprint
  z = r - sigma
  hsigma = (1 - 0.5*eta)/(1-eta)**3 - 1
  density = 3/4/pi*eta
  rhs = (1-eta)**4/(1+4*eta+4*eta**2-4*eta**3+eta**4)

  a0 = xnew[0]
  a1 = 0 # this parameter is constrained via integral
  a2 = xnew[1]
  a3 = xnew[2]
  a4 = 0
  a5 = xnew[3]

  int_h0 = 4*pi*hsigma*(2 + sigma*a0*(2 + sigma*a0))/a0**3 - 4/3*pi*sigma**3

  int_h1_over_a1 = 4*pi*a2* \
      (sigma**2*a2**4 + 2*a2**2* \
         (-1 + sigma*a3*(2 + sigma*a3)) + a3**2*(6 + sigma*a3*(4 + sigma*a3))) \
         /(a2**2 + a3**2)**3

  int_h2_over_a4 = 4*pi*a5*(6 + sigma*a5*(4 + sigma*a5))/a5**4

  A = (rhs-1)/density - int_h0
  B = int_h2_over_a4
  C = int_h1_over_a1

  a1 = (A - B*hsigma*a0 + B*hsigma + B*hsigma**2)/(C - B*a2)

  a4 = hsigma*a0 - a1*a2 - hsigma - hsigma**2

  f0 = hsigma
  j0 = exp(-a0*z)

  f1 = a1
  j1 = sin(a2*z)*exp(-a3*z)

  f2 = a4
  j2 = z*exp(-a5*z)

  return 1 + f0*j0 + f1*j1 + f2*j2

def evalg_cubic(xnew, eta, r):
  global toprint
  z = r - sigma
  hsigma = (1 - 0.5*eta)/(1-eta)**3 - 1
  density = 3/4/pi*eta
  rhs = (1-eta)**4/(1+4*eta+4*eta**2-4*eta**3+eta**4)

  a0 = xnew[0]
  a1 = 0
  a2 = xnew[1]
  a3 = 0
  a4 = xnew[2]
  a5 = xnew[3]
  a6 = xnew[4]

  # a6 = a0
  # a4 = a0
  # a2 = a0

  int_h0 = 4*pi*hsigma*(2 + sigma*a0*(2 + sigma*a0))/a0**3 - 4/3*pi*sigma**3

  int_h1_over_a1 = 4*pi*(6 + sigma*a2*(4 + sigma*a2))

  int_h2_over_a3 = 8*pi*(12 + sigma*a4*(6 + sigma*a4))/a4**5

  int_h3 = a5*24*pi*(20 + sigma*a6*(8 + sigma*a6))/a6**6

  A = (rhs-1)/density - int_h0 - int_h3
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

  f3 = a5
  j3 = z**3*exp(-a6*z)

  return 1 + f0*j0 + f1*j1 + f2*j2 + f3*j3

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

chi2 = sum(dist2(x)**2)
print "beginning least squares fit, chi^2 initial: %g" %chi2
vals, mesg = leastsq(dist2, x)
chi2 = sum(dist2(vals)**2)
print "original fit complete, chi^2: %.3f" % chi2

toprint = True
for i in range(len(x)):
  print "vals[%i]: %.3f\t x[%i]: %g" %(i, vals[i], i, x[i])

g = dist(vals)
gdifference = dist2(vals)

# for i in g:
#     print dist(vals, ind)[i]

for i in range(len(ff)):
  figure(1)
  plot(r_mc, g[i*len(r):(i+1)*len(r)], colors[i]+'--',label='g at filling fraction %.2f'%ff[i])
  hsigma = (1 - 0.5*ff[i])/(1-ff[i])**3 - 1
  density = 4/3*pi*ff[i]
  rhs = (1-ff[i])**4/(1+4*ff[i]+4*ff[i]**2-4*ff[i]**3+ff[i]**4)/3
  #integral = hsigma*(1/a + x[0]*x[1]/())



  #print density, integral, rhs
  #print "ff: %.2f\t thing: %g" %(ff[i], 1 - rho*integral - rhs)
  figure(2)
  #plot(r_mc, gdifference[i*len(r):(i+1)*len(r)], colors[i]+'--')
  plot(r_mc, g[i*len(r):(i+1)*len(r)] - ghs[i], colors[i]+'-')
  # calculating integral:
  #mc:
  r_mc, ghs[i]
  integrand_mc = 4*pi*r_mc*r_mc*ghs[i]
  integrand_ours = 4*pi*r_mc*r_mc*g[i*len(r):(i+1)*len(r)]
  integral_mc = sum(integrand_mc)/len(integrand_mc)*(r_mc[2]-r_mc[1]) - 4/3*pi*sigma**3
  integral_ours = sum(integrand_ours)/len(integrand_ours)*(r_mc[2]-r_mc[1]) - 4/3*pi*sigma**3
  print("Int_mc: %6.3f, Int_ours: %6.3f, Diff: %6.3f" %(integral_mc, integral_ours, integral_ours-integral_mc))



  #plot(r_mc, numpy.abs(numpy.asarray(ghsconcatenated[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'-')



figure(1)
xlim(2,6.5)
ylim(0., 3.5)
xlabel(r"$r/R$")
ylabel("$g(r)$")
legend(loc='best').get_frame().set_alpha(0.5)

savefig("figs/ghs-g2.pdf")


figure(2)
xlim(2,6.5)
ylim(-.25, .25)
xlabel(r"$r/R$")
ylabel("|ghs - g|")
#legend(loc='best').get_frame().set_alpha(0.5)
savefig("figs/ghs-g-ghs2.pdf")

# figure(3)
# for eta in [.5, .6, .7, .8]:
#   gsigma = (1-eta/2)/(1-eta)**3
#   plot(r_mc, evalg(vals, eta, r), label='eta %g  gsig %g'%(eta, gsigma))
axhline(y=0)
xlim(2,6.5)
legend(loc='best')
show()
