#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
#matplotlib.use('Agg')

import pylab, numpy, sys
import Scientific.Functions.LeastSquares as ls
#if len(sys.argv) != 2:
#    print("Usage:  " + sys.argv[0] + " out-filename.pdf")
#    exit(1)

#pdffilename = sys.argv[1]

pylab.figure(1)
pylab.title('$g_{HS}(r)$')
pylab.axvline(x=1, color='k', linestyle=':')
pylab.axhline(y=1, color='k', linestyle=':')

#pylab.figure(2)
#pylab.axvline(x=1, color='k', linestyle=':')
#pylab.axhline(y=0, color='k', linestyle='-')

toFit = 4
numParams = 8
x = [1]*numParams
x[0] = 1.15
x[1] = 6.8

x[2] = .5
x[3] = 6.56
x[4] = 2.1

x[5] = -.4
x[6] = 4.8
x[7] = 3.02

def dist(x, ind):
    gsig, r = ind
    gsig -= 1
    # function with x[i] as constants to be determined
    h0 = x[0]*gsig
    f0 = numpy.exp(-x[1]*r)

    h1 = x[2]*gsig
    f1 = numpy.sin(x[3]*r) * numpy.exp(-x[4]*r)
    
    h2 = x[5]*gsig**(2)
    f2 = numpy.sin(x[6]*r) * numpy.exp(-x[7]*r)
   
    g = 1 + h0*f0 + h1*f1 + h2*f2
    return g

def distWeighted(x, ind):
    gsig, r = ind
    gsig -= 1
    # function with x[i] as constants to be determined
    h0 = x[0]*gsig
    f0 = numpy.exp(-x[1]*r)

    h1 = x[2]*gsig
    f1 = numpy.sin(x[3]*r) * numpy.exp(-x[4]*r)
    
    h2 = x[5]*gsig**(2)
    f2 = numpy.sin(x[6]*r) * numpy.exp(-x[7]*r)
   
    g = 1 + h0*f0 + h1*f1 + h2*f2
    return g*r**(-4)
  
def plotdatag(x, ind):
    gsig, r = ind
    #gsig -= 1

    # function with x[i] as constants to be determined
    h0 = x[0]*gsig**2
    f0 = numpy.exp(-x[1]*r - x[2])

    h1 = x[3]*gsig
    f1 = numpy.cos(x[4]*r - x[5]) * numpy.exp(-x[6]*r)
    
    h2 = x[7]*gsig**(-1)
    f2 = numpy.cos(x[8]*r - x[9]) * numpy.exp(-x[10]*r)

    h3 = x[11]*gsig**(2)
    f3 = numpy.cos(x[12]*r - x[13]) * numpy.exp(-x[14]*r)
    
#    h4 = x[15]*gsig**(-2)
#    f4 = numpy.cos(x[16]*r - x[17]) * numpy.exp(-x[18]*r)
       
    g = 1 + h0*f0 + h1*f1 + h2*f2 + h3*f3 #+ h4*f4
    return g, h0*f0, h1*f1, h2*f2, h3*f3
    
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

colors = ['r', 'g', 'b', 'c']
ff = [.1008, .2, .3, .4]

# READ DATA
ghs = [0]*4
gsig = [0]*4
i = 0
while (i < 4):
    r_mc, ghs[i] = read_ghs(ff[i])

    pylab.figure(1)
    pylab.plot(r_mc, ghs[i], colors[i]+"-",label='filling fraction %.1f'%ff[i])

    gsig[i] = ghs[i].max()
    r = (r_mc)/2 - 1 # "diameter units", shifted to x = 0 at d = 1
    i += 1

g2 = numpy.concatenate((ghs[0], ghs[1], ghs[2], ghs[3])) # make g2 = g(sig, r) = [g(sig1, r1), g(sig1, r2), ..., g(sig2, r2), ...]

ind = [0]*len(r)*len(gsig)
j = 0
while (j < len(gsig)):   # makes ind an array of pairs, ind = [(sig1, r1), (sig1, r2), ..., (sig2, r2), ...]
    i = 0
    while (i < len(r)):
        ind[i + j*len(r)] = (gsig[j], r[i])
        i += 1
    j += 1

vals = [0]*numParams
valsw = [0]*numParams
#vals, cov = curve_fit(dist, r, ghs, x, sig)

i = 0
data = [0]*len(r)*len(gsig)
dataW = [0]*len(r)*len(gsig)

while (i < len(r)*len(gsig)):
    data[i] = (ind[i], g2[i])
    dataW[i] = (ind[i], g2[i], (ind[i][1]+1)**4)
#    dataW[i] = ((ind[i][0], ind[i][1]+1), g2[i]*(ind[i][1]+1)**(-4))
#    print "ind[i] + (0,1): " + str((ind[i][0], ind[i][1]+1)) + "r: " + str((ind[i][1]+1)) + ", r^-4: " + str((ind[i][1]+1)**(-4)) + ", gW: " + str(g2[i]*(ind[i][1]+1)**(-4))
    i += 1

vals, cov = ls.leastSquaresFit(dist, x, data)
print "original fit complete, cov: " + str(cov)
i = 0
while (i < len(vals)):
    print 'x[' + str(i) + ']: ' + str(vals[i])
    i += 1

valsw, cov = ls.leastSquaresFit(dist, x, dataW)
print "weighted fit complete, cov: " + str(cov)
i = 0
while (i < len(vals)):
    print 'x[' + str(i) + ']: ' + str(valsw[i])
    i += 1


i = 0
g = numpy.zeros(len(r)*len(gsig))
gw = numpy.zeros(len(r)*len(gsig))
while (i < len(r)*len(gsig)):
    #g[i], g0[i], g1[i], g2[i], g3[i] = plotdatag(vals, ind[i])
    g[i] = dist(vals, ind[i])
    gw[i] = dist(valsw, ind[i])
#    gw[i] = distWeighted(valsw, (ind[i][0], ind[i][1]+1))*(ind[i][1] + 1)**(4)
    i += 1

for (color,i) in [('r', 0), ('g', 1), ('b', 2), ('c', 3)]:
#for (color,i) in [('r', 0), ('g', 1), ('b', 2)]:
    pylab.figure(1)
    pylab.plot(r_mc, g[i*len(r):(i+1)*len(r)], color+'--')
#    pylab.plot(r_mc, gw[i*len(r):(i+1)*len(r)], color+'--.')

#    pylab.plot(r_mc, g0[i*len(r):(i+1)*len(r)]+g3[i*len(r):(i+1)*len(r)]+1, color+'--.')
    pylab.figure(2)
    pylab.plot(r_mc, numpy.abs(numpy.asarray(g[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'-')
#    pylab.plot(r_mc, numpy.abs(numpy.asarray(gw[i*len(r):(i+1)*len(r)]) - ghs[i]), color+'--')


pylab.figure(1)
pylab.xlim(2,6.5)
#pylab.ylim(0.5, 3.5)
pylab.xlabel(r"$r/R$")
pylab.ylabel("$g(r)$")
pylab.legend(loc='best').get_frame().set_alpha(0.5)
#pylab.savefig(pdffilename)



pylab.figure(2)
pylab.xlim(2,6.5)
pylab.ylim(0, .1)
pylab.xlabel(r"$r/R$")
pylab.ylabel("|ghs - g|")
#pylab.legend(loc='best').get_frame().set_alpha(0.5)

pylab.show()
