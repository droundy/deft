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

mcoffset = 13
n0 = dftdata[:,6]
nA = dftdata[:,8]
nAmc = mcdata[:,11]
n0mc = mcdata[:,10]

pylab.figure(figsize=(8, 12))
pylab.subplots_adjust(hspace=0.001)

n_plt = pylab.subplot(3,1,3)
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
    pylab.ylim(0.1,0.8)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
    pylab.ylim(0.17,0.43)
if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
    pylab.ylim(0.04,0.19)
pylab.xlabel("position")
#pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

#pylab.twinx()

stop_here = int(dft_len - 1/dft_dr)
print stop_here
start_here = int(2.5/dft_dr)
off = 1
me = 30

A_plt = pylab.subplot(3,1,1)
A_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,2+2*off]/nAmc,"b-",label="$g_\sigma^A$ MC")
A_plt.plot(dftdata[start_here:stop_here,0],dftdata[start_here:stop_here,5],"g--",markevery=me,label="$g_\sigma^A$ (White Bear)")
A_plt.plot(wbm2data[start_here:stop_here,0],wbm2data[start_here:stop_here,5],"m--",markevery=me,label="$g_\sigma^A$ (mark II)",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.plot(dftdata[:,0],dftdata[:,7],"rx--",markevery=me,label="Gross",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
A_plt.legend(loc=1, ncol=1).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
#matplotlib.ticks(arange(0.5, 1.5, 3.5))
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
#     pylab.ylim(2.7,4.8)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
#     pylab.ylim(1.0,4.20)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
#     pylab.ylim(0.95,1.70)

S_plt = pylab.subplot(3,1,2)
S_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,3+2*off]/n0mc,"b-",label="$g_\sigma^S$ MC")
S_plt.plot(dftdata[start_here:stop_here,0],dftdata[start_here:stop_here,3],
           "g--",markevery=me,label="$g_\sigma^S$ (White Bear)")
S_plt.plot(wbm2data[start_here:stop_here,0],wbm2data[start_here:stop_here,3],"m--",markevery=me,label="$g_\sigma^S$ (mark II)",
           markerfacecolor='none',markeredgecolor='green', markeredgewidth=1)
S_plt.plot(dftdata[:,0],dftdata[:,4],"rx--",markevery=me,label="Yu and Wu")
S_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)
pylab.xlim(0,12)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
#     pylab.ylim(2.3,5.50)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
#     pylab.ylim(1.50,3.60)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
#     pylab.ylim(1.10,1.60)

pylab.ylabel("filling fraction")

xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
pylab.setp(xticklabels, visible=False)

pylab.savefig(sys.argv[5])

pylab.show()

# pylab.figure(figsize=(8,8))
# pylab.subplots_adjust(hspace=0.001)

# Sgnn_plt = pylab.subplot(2,1,1)
# Sgnn_plt.plot(dftdata[start_here:stop_here,0],dftdata[start_here:stop_here,1]*nA[start_here:stop_here]*(4*numpy.pi/3)**2*dftdata[start_here:stop_here,5],
#            "g+--",markevery=me,label="$nn_Ag_\sigma^A$ (White Bear)")
# Sgnn_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,2+2*off]*mcdata[:,1]*(4*numpy.pi/3)**2,"g-",label="$nn_Ag_\sigma^A$ MC")
# Sgnn_plt.plot(dftdata[start_here:stop_here,0],(dftdata[start_here:stop_here,1]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,7],
#            "rx--",markevery=me,label="$nn_Ag_\sigma^A$ (Gross)")
# Sgnn_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)
# #It seems like with the Gross here we do n*n*g and not n*nA*g, why?

# Agnn_plt = pylab.subplot(2,1,2)
# Agnn_plt.plot(dftdata[start_here:stop_here,0],(n0[start_here:stop_here]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,3],
#            "g+--",markevery=me,label="$n_0^3g_\sigma^S$ (White Bear)")
# Agnn_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,3+2*off]*n0mc*(4*numpy.pi/3)**2,"g-",label="$n_0^2g_\sigma^S$ MC")
# Agnn_plt.plot(dftdata[start_here:stop_here,0],(n0[start_here:stop_here]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,4],
#            "rx--",markevery=me,label="$n_0^2g_\sigma^S$ (YuWu)")
# Agnn_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)

# pylab.xlim(0,12)
# xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
# pylab.setp(xticklabels, visible=False)


# pylab.savefig("gnn-" + sys.argv[5])


# div = len(dftdata[:,0])
# sep = dftdata[1,0] - dftdata[0,0]
# L = sep*div
# mid = int (div/2.0)
# midsec = int (div/2.0) - 20
# Nb = 0
# for i in dftdata[midsec:mid,0]
#   Ni = dftdata[i,0]*sep
#   Nb = Nb +Ni

# nb = Nb/(20*20*20*sep*sep*sep) 

# Ntot = 0
# for j in dftdata[:,0]:
#   nj = dftdata[j,0]*sep
#   Ntot = Ntot + nj

# Ntot = Ntot/(sep*div*sep*div)

# Nextra = Ntot - L*nb
#want to put in terms of per Area and also put the whole thing in th cpp file.  The abovefor j in dftdata[:,0] is 
#quesionable as to what that will go through


# N = 0
# Totvol = 0
# for i in dftdata[:,0]:
#   vol = 4*numpy.pi/3*dftdata[i+1,0]*dftdata[i+1,0]*dftdata[i+1,0] - 4*numpy.pi/3*dftdata[i,0]*dftdata[i,0]*dftdata[i,0]
#   N = N + dftdata[i,1]*vol
#   Totvol = Totvol +vol

# NperVol = N/Totvol

# #print("For " + sys.argv[2] + " total number per volume = %f. Multiply this by your total volume to get the same filling fraction.",(N/Totvol),)jj 

# strFile = "figs/FillingFracInfo.dat"
# file = open(strFile , "a")
# file.write("For " + sys.argv[2] + " total number per volume = " + str(NperVol) + ".\n")
# file.close()
