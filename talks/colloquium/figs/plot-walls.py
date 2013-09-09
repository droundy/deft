#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

xmax = 2.5
xmin = -0.4

def plotit(dftdata, mcdata):
    dft_len = len(dftdata[:,0])
    dft_dr = dftdata[2,0] - dftdata[1,0]

    mcdata = numpy.insert(mcdata,0,0,0)
    mcdata[0,0]=-10

    mcoffset = 10/2
    offset = -3/2
    n0 = dftdata[:,6]
    nA = dftdata[:,8]
    nAmc = mcdata[:,11]
    n0mc = mcdata[:,10]

    pylab.figure(figsize=(6, 6))
    pylab.subplots_adjust(hspace=0.001)

    n_plt = pylab.subplot(3,1,3)
    n_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,1]*4*numpy.pi/3,"b-",label='$n$ Monte Carlo')
    n_plt.plot(dftdata[:,0]/2+offset,dftdata[:,1]*4*numpy.pi/3,"b--",label='$n$ DFT')
    n_plt.legend(loc='best', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
    n_plt.yaxis.set_major_locator(pylab.MaxNLocator(6,steps=[1,5,10],prune='upper'))
    pylab.ylim(ymin=0)
    pylab.xlim(xmin, xmax)
    pylab.xlabel("$z/\sigma$")
    pylab.ylabel("$n(\mathbf{r})$")
    n_plt.axvline(x=0, color='k', linestyle=':')


    n = len(mcdata[:,0])

#pylab.twinx()

    stop_here = int(dft_len - 1/dft_dr)
    print stop_here
    start_here = int(2.5/dft_dr)
    off = 1
    me = 55

    A_plt = pylab.subplot(3,1,1)
    A_plt.axvline(x=0, color='k', linestyle=':')
    A_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,2+2*off]/nAmc,"r-",label="$g_\sigma^A$ Monte Carlo")
    A_plt.plot(dftdata[start_here:stop_here,0]/2+offset,dftdata[start_here:stop_here,5],"ro",markevery=me,label="$g_\sigma^A$ this work")
    A_plt.plot(dftdata[:,0]/2+offset,dftdata[:,7],"rx",markevery=me,label="Gross",
               markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
    A_plt.legend(loc='best', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
    A_plt.yaxis.set_major_locator(pylab.MaxNLocator(integer=True,prune='upper'))
    pylab.ylim(ymin=0)
    pylab.ylabel("$g_\sigma^A$")
    pylab.xlim(xmin, xmax)

    n0mc[0]=1
    mcdata[0,10]=1
    S_plt = pylab.subplot(3,1,2)
    S_plt.axvline(x=0, color='k', linestyle=':')
    S_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,3+2*off]/n0mc,"g-",label="$g_\sigma^S$ Monte Carlo")
    S_plt.plot(dftdata[:,0]/2+offset,dftdata[:,4],"gx",markevery=me/2,label="Yu and Wu")
    S_plt.legend(loc='best', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
    #pylab.ylim(ymax=12)
    S_plt.yaxis.set_major_locator(pylab.MaxNLocator(5,integer=True,prune='upper'))
    pylab.xlim(xmin, xmax)
    pylab.ylim(ymin=0)
    pylab.ylabel("$g_\sigma^S$")

    xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
    pylab.setp(xticklabels, visible=False)

mcdata10 = numpy.loadtxt('../../papers/contact/figs/mc-walls-20-196.dat')
dftdata10 = numpy.loadtxt('../../papers/contact/figs/wallsWB-0.10.dat')

mcdata40 = numpy.loadtxt('../../papers/contact/figs/mc-walls-20-817.dat')
dftdata40 = numpy.loadtxt('../../papers/contact/figs/wallsWB-0.40.dat')

plotit(dftdata10, mcdata10)
pylab.savefig('figs/walls-10.pdf', transparent=True)

plotit(dftdata40, mcdata40)
pylab.savefig('figs/walls-40.pdf', transparent=True)
