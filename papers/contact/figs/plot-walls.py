#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 6:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat wbm2-filename.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])
wbm2data = numpy.loadtxt(sys.argv[4])

mcoffset = 8.5

pylab.plot(mcdata[:,0]+mcoffset,mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0]-3,dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT density')
pylab.plot(wbtdata[:,0]-3,wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT density')
pylab.plot(wbm2data[:,0]-3,wbm2data[:,1]*4*numpy.pi/3,"c--",label='WB mark II density')

me = 20

pylab.plot(dftdata[:,0]-3,dftdata[:,6]*4*numpy.pi/3,"c--",markevery=me,label="$n_0$")

pylab.plot(mcdata[:,0]+mcoffset,mcdata[:,4]*4*numpy.pi/3,"r-",label="MC $n^A$")
pylab.plot(mcdata[:,0]+mcoffset,mcdata[:,5]*4*numpy.pi/3,"g-",label="MC $n^S$")

pylab.plot(dftdata[:,0]-3,dftdata[:,3]*4*numpy.pi/3,"g+--",markevery=me,label="$n^{S}_{contact}$")
pylab.plot(dftdata[:,0]-3,dftdata[:,4]*4*numpy.pi/3,"gx--",markevery=me,label="Yu and Wu")
pylab.plot(dftdata[:,0]-3,dftdata[:,5]*4*numpy.pi/3,"ro--",markevery=me,label="DFT at sphere")
pylab.plot(wbm2data[:,0]-3,wbm2data[:,3]*4*numpy.pi/3,"go--",markevery=me,label="$n^{S}_{contact}$ (mark II)",
           markerfacecolor='none',markeredgecolor='green', markeredgewidth=1)
pylab.plot(wbm2data[:,0]-3,wbm2data[:,5]*4*numpy.pi/3,"r+--",markevery=me,label="DFT at sphere (mark II)",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
pylab.plot(dftdata[:,0]-3,dftdata[:,7]*4*numpy.pi/3,"rx--",markevery=me,label="Gross",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)

pylab.xlabel("z")
pylab.ylabel("filling fraction")
pylab.legend(loc='lower right', ncol=2).get_frame().set_alpha(0.5)
pylab.xlim(-1,8)
pylab.ylim(ymin=0, ymax=1.2*dftdata[:,1].max()*4*numpy.pi/3)

pylab.savefig(sys.argv[5])

NperVol = 0
N = 0
Totvol = 0
for i in range(len(dftdata[:,0])-1):
    if dftdata[i,1] > 0.001:
        vol = 4*numpy.pi/3*dftdata[i+1,0]*dftdata[i+1,0]*dftdata[i+1,0] - 4*numpy.pi/3*dftdata[i,0]*dftdata[i,0]*dftdata[i,0]
        N = N + dftdata[i,1]*vol
        Totvol = Totvol +vol
    
NperVol = N/Totvol

#print("For " + sys.argv[2] + " total number per volume = %f. Multiply this by your total volume to get the same filling fraction.",(N/Totvol),) 

strFile = "figs/FillingFracInfo.dat"
file = open(strFile , "a")
file.write("For " + sys.argv[2] + " total number per volume = " +  str(NperVol) +".\n")
file.close()
