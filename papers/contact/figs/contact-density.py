#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
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
n0 = dftdata[:,6]
nA = dftdata[:,8]
nAmc = mcdata[:,11]
n0mc = mcdata[:,10]

pylab.figure(figsize=(8, 8))
pylab.subplots_adjust(hspace=0.001)

n_plt = pylab.subplot(3,1,3)
n_plt.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC $n$')
n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT $n$')
n_plt.plot(wbtdata[:,0],wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT $n$')
n_plt.plot(wbm2data[:,0],wbm2data[:,1]*4*numpy.pi/3,"c--",label='Mark II $n$')
n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")
n_plt.plot(mcdata[:,0],n0mc*4*numpy.pi/3,"c-",label="MC $n_0$")
n_plt.plot(mcdata[:,0],nAmc*4*numpy.pi/3,"m-",label="MC $n_A$")
n_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=2).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
n = len(mcdata[:,0])
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
    pylab.ylim(0.1,0.8)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
    pylab.ylim(0.17,0.43)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
    pylab.ylim(0.04,0.19)
pylab.xlabel("position")
pylab.ylabel("filling fraction")

#pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

#pylab.twinx()

off = 2
me = 3
A_plt = pylab.subplot(3,1,1)
A_plt.plot(mcdata[:,0],mcdata[:,2+2*off]/nAmc,"r-",label="MC $g_\sigma^A$")
A_plt.plot(dftdata[:,0],dftdata[:,5],"ro--",label="$g_\sigma^A$ (White Bear)")
A_plt.plot(wbm2data[:,0],wbm2data[:,5],"r+--",markevery=me,label="$g_\sigma^A$ (mark II)",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.plot(dftdata[:,0],dftdata[:,7],"rx--",markevery=me,label="Gross",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.legend(loc=1, ncol=1).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
#matplotlib.ticks(arange(0.5, 1.5, 3.5))
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
    pylab.ylim(2.7,4.8)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
    pylab.ylim(1.0,4.20)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
    pylab.ylim(0.95,1.70)

sphere_end = int(dft_len - 1/dft_dr)
print sphere_end

S_plt = pylab.subplot(3,1,2)
S_plt.plot(mcdata[:,0],mcdata[:,3+2*off]/n0mc,"g-",label="MC $g_\sigma^S$")
S_plt.plot(dftdata[0:sphere_end,0],dftdata[0:sphere_end,3],
           "g+--",label="$g_\sigma^S$ (White Bear)")
S_plt.plot(wbm2data[:,0],wbm2data[:,3],"go--",markevery=me,label="$g_\sigma^S$ (mark II)",
           markerfacecolor='none',markeredgecolor='green', markeredgewidth=1)
S_plt.plot(dftdata[:,0],dftdata[:,4],"gx--",label="Yu and Wu")
S_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
    pylab.ylim(2.3,5.50)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
    pylab.ylim(1.50,3.60)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
    pylab.ylim(1.10,1.60)

pylab.ylabel("filling fraction")

xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
pylab.setp(xticklabels, visible=False)

pylab.savefig(sys.argv[5])

pylab.show()



N = 0
Totvol = 0
for i in dftdata[:,0]:
  vol = 4*numpy.pi/3*dftdata[i+1,0]*dftdata[i+1,0]*dftdata[i+1,0] - 4*numpy.pi/3*dftdata[i,0]*dftdata[i,0]*dftdata[i,0]
  N = N + dftdata[i,1]*vol
  Totvol = Totvol +vol

NperVol = N/Totvol

#print("For " + sys.argv[2] + " total number per volume = %f. Multiply this by your total volume to get the same filling fraction.",(N/Totvol),)jj 

strFile = "figs/FillingFracInfo.dat"
file = open(strFile , "a")
file.write("For " + sys.argv[2] + " total number per volume = " + str(NperVol) + ".\n")
file.close()

