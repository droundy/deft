#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 6:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat wb-m2.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])
wbm2data = numpy.loadtxt(sys.argv[4])

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

pylab.figure(figsize=(7, 6))
pylab.subplots_adjust(hspace=0.001)

n_plt = pylab.subplot(3,1,3)
#n_plt.axvline(x=3, color='k', linestyle=':')
n_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,1]*4*numpy.pi/3,"b-",label='$n$ MC')
n_plt.plot(dftdata[:,0]/2+offset,dftdata[:,1]*4*numpy.pi/3,"b--",label='$n$ DFT')
#n_plt.plot(wbtdata[:,0]/2+offset,wbtdata[:,1]*4*numpy.pi/3,"c-.",label='$n$ WBT')
#n_plt.plot(wbm2data[:,0]/2+offset,wbm2data[:,1]*4*numpy.pi/3,"m--",label='$n$ Mark II')
#n_plt.plot(dftdata[:,0]/2+offset,n0*4*numpy.pi/3,"g--",label="$n_0$")
#n_plt.plot(dftdata[:,0]/2+offset,nA*4*numpy.pi/3,"g-.",label="$n_A$")
#n_plt.plot(dftdata[:,0]/2+offset,n0*4*numpy.pi/3,"c-.",label="$n_0$")
#n_plt.plot(dftdata[:,0]/2+offset,nA*4*numpy.pi/3,"m--",label="$n_A$")
if (dftdata[dft_len-2,1]*4*numpy.pi/3 < 0.15 and dftdata[dft_len-2,1]*4*numpy.pi/3 > 0.05):
    n_plt.legend(loc='lower right', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
else:
    n_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
n_plt.yaxis.set_major_locator(pylab.MaxNLocator(6,steps=[1,5,10],prune='upper'))
pylab.ylim(ymin=0)
pylab.xlim(-1/2,9/2)
pylab.xlabel("$r$/$\sigma$")
pylab.ylabel("filling fraction")
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
A_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,2+2*off]/nAmc,"r-",label="$g_\sigma^A$ MC")
A_plt.plot(dftdata[start_here:stop_here,0]/2+offset,dftdata[start_here:stop_here,5],"ro-",markevery=me,label="$g_\sigma^A$ this work")
#A_plt.plot(wbm2data[start_here:stop_here,0]/2+offset,wbm2data[start_here:stop_here,5],"rx-",markevery=me,label="$g_\sigma^A$ (mark II)",
#           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.plot(dftdata[:,0]/2+offset,dftdata[:,7],"rx--",markevery=me,label="Gross",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.legend(loc='lower right', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
A_plt.yaxis.set_major_locator(pylab.MaxNLocator(integer=True,prune='upper'))
pylab.ylim(ymin=0)
#pylab.ylim(ymax=(dftdata[20,5]+1))
pylab.ylabel("$g^A$")
pylab.xlim(-1/2,9/2)

n0mc[0]=1
mcdata[0,10]=1
S_plt = pylab.subplot(3,1,2)
S_plt.axvline(x=0, color='k', linestyle=':')
S_plt.plot(mcdata[:,0]/2+mcoffset,mcdata[:,3+2*off]/n0mc,"g-",label="$g_\sigma^S$ MC")
S_plt.plot(dftdata[start_here:stop_here,0]/2+offset,dftdata[start_here:stop_here,3],
           "go--",markevery=me,label="$g_\sigma^S$ this work")
#S_plt.plot(wbm2data[start_here:stop_here,0]/2+offset,wbm2data[start_here:stop_here,3],"g--",markevery=me,label="$g_\sigma^S$ (mark II)",
#           markerfacecolor='none',markeredgecolor='green', markeredgewidth=1)
S_plt.plot(dftdata[:,0]/2+offset,dftdata[:,4],"gx--",markevery=me,label="Yu and Wu")
if (dftdata[dft_len-2,1]*4*numpy.pi/3 < 0.15 and dftdata[dft_len-2,1]*4*numpy.pi/3 > 0.05):
    S_plt.legend(loc='lower right', ncol=1).draw_frame(False) #.get_frame().set_alpha(0.5)
else:
    S_plt.legend(loc='upper right', ncol=2).draw_frame(False) #get_frame().set_alpha(0.5)
    #pylab.ylim(ymax=12)
S_plt.yaxis.set_major_locator(pylab.MaxNLocator(5,integer=True,prune='upper'))
pylab.xlim(-1/2,9/2)
pylab.ylim(ymin=0)
pylab.ylabel("$g^S$")

xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
pylab.setp(xticklabels, visible=False)

pylab.savefig(sys.argv[5])

pylab.show()
