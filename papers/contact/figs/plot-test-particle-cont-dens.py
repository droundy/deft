#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, math, scipy

from scipy.interpolate import interp1d

if len(sys.argv) != 7:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat wb-m2.dat test-particle-filename.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])
wbm2data = numpy.loadtxt(sys.argv[4])
testdata = numpy.loadtxt(sys.argv[5])

dft_len = len(dftdata[:,0])
dft_dr = dftdata[2,0] - dftdata[1,0]
test_len = len(testdata[:,0])

# for x in a[:]: # make a slice copy of the entire list
# ...    if len(x) > 6: a.insert(0, x)
# ...
# >>> a
# ['defenestrate', 'cat', 'window', 'defenestrate']

pos = 0
minpos = [0]*test_len
minposmc = [0]*test_len
new_gA = [0]*test_len
new_gmc = [0]*test_len
new_nAmc = [0]*test_len
new_nA = [0]*test_len
mcoffset = 13

for i in range(test_len):
    new_nA[i] = numpy.interp(testdata[i,0], dftdata[:,0], dftdata[:,8])
    new_gA[i] = numpy.interp(testdata[i,0], dftdata[:,0], dftdata[:,5])
    new_nAmc[i] = numpy.interp(testdata[i,0], mcdata[:,0] + mcoffset, mcdata[:,11])
    new_gmc[i] = numpy.interp(testdata[i,0], mcdata[:,0] + mcoffset, mcdata[:,4])/new_nAmc[i]


n0 = dftdata[:,6]
nA = dftdata[:,8]
nAmc = mcdata[:,11]
n0mc = mcdata[:,10]

pylab.figure(figsize=(8, 8))
pylab.subplots_adjust(hspace=0.001)

n_plt = pylab.subplot(2,1,2)
n_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,1]*4*numpy.pi/3,"b-",label='$n$ MC')
n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"g--",label='$n$ DFT')
n_plt.plot(wbtdata[:,0],wbtdata[:,1]*4*numpy.pi/3,"c-.",label='$n$ WBT')
n_plt.plot(wbm2data[:,0],wbm2data[:,1]*4*numpy.pi/3,"m--",label='$n$ Mark II')
n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"g--",label="$n_0$")
n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"g-.",label="$n_A$")
n_plt.plot(mcdata[:,0]+mcoffset,n0mc*4*numpy.pi/3,"bx",label="$n_0$ MC")
n_plt.plot(mcdata[:,0]+mcoffset,nAmc*4*numpy.pi/3,"b+",label="$n_A$ MC")
n_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=2).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
n = len(mcdata[:,0])
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
    pylab.ylim(0.1,2.8)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
    pylab.ylim(0.17,1.19)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
    pylab.ylim(0.04,0.16)
    n_plt.yaxis.set_ticks(pylab.arange(0.0, 0.15, 0.02))
    n_plt.legend(loc='lower right', ncol=2).get_frame().set_alpha(0.5)
pylab.ylim(ymin=0)
pylab.xlabel("position")
#pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

#pylab.twinx()

stop_here = int(dft_len - 1/dft_dr)
print stop_here
start_here = int(2.5/dft_dr)
off = 1
me = 30

T_plt = pylab.subplot(2,1,1)
T_plt.plot(testdata[:,0],testdata[:,1]/new_nA[:],"bo",label="$g_\sigma^A$ Test-particle")
T_plt.plot(testdata[:,0],new_gA[:],"ro",label="$g_\sigma^A$ DFT")
T_plt.plot(testdata[:,0],new_gmc[:],"go",label="$g_\sigma^A$ MC")
T_plt.legend(loc='lower right', ncol=1).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
#pylab.ylim(-5,20)

pylab.ylabel("filling fraction")

xticklabels = T_plt.get_xticklabels()
pylab.setp(xticklabels, visible=False)

pylab.savefig(sys.argv[6])

pylab.show()
