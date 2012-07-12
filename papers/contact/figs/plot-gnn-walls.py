#!/usr/bin/python

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

stop_here = int(dft_len - 1/dft_dr)
print stop_here
start_here = int(2.5/dft_dr)
off = 0
me = 30
n = len(mcdata[:,0])

pylab.figure(figsize=(8,6))
pylab.subplots_adjust(hspace=0.001)

Agnn_plt = pylab.subplot(2,1,1)
Agnn_plt.plot(dftdata[start_here:stop_here,0],dftdata[start_here:stop_here,1]*nA[start_here:stop_here]*(4*numpy.pi/3)**2*dftdata[start_here:stop_here,5],
           "g+--",markevery=me,label="$nn_Ag_\sigma^A$ (White Bear)")
Agnn_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,2+2*off]*mcdata[:,1]*(4*numpy.pi/3)**2,"g-",label="$nn_Ag_\sigma^A$ MC")
Agnn_plt.plot(dftdata[start_here:stop_here,0],(dftdata[start_here:stop_here,1]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,7],
           "rx--",markevery=me,label="$nn_Ag_\sigma^A$ (Gross)")
Agnn_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)
#It seems like with the Gross here we do n*n*g and not n*nA*g, why?
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
#     pylab.ylim(0.3,0.980) 
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
#     pylab.ylim(0.090,0.360)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
#     pylab.ylim(0.000,0.034)

Sgnn_plt = pylab.subplot(2,1,2)
Sgnn_plt.plot(dftdata[start_here:stop_here,0],(n0[start_here:stop_here]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,3],
           "g+--",markevery=me,label="$n_0^2g_\sigma^S$ (White Bear)")
Sgnn_plt.plot(mcdata[:,0]+mcoffset,mcdata[:,3+2*off]*n0mc*(4*numpy.pi/3)**2,"g-",label="$n_0^2g_\sigma^S$ MC")
Sgnn_plt.plot(dftdata[start_here:stop_here,0],(n0[start_here:stop_here]*4*numpy.pi/3)**2*dftdata[start_here:stop_here,4],
           "rx--",markevery=me,label="$n_0^2g_\sigma^S$ (YuWu)")
Sgnn_plt.legend(loc=1, bbox_to_anchor=[1.0, 1.0], ncol=1).get_frame().set_alpha(0.5)

pylab.xlim(0,12)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
#     pylab.ylim(0.270,0.980)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
#     pylab.ylim(0.120,0.380)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
#     pylab.ylim(0.000,0.030)
#xticklabels = A_plt.get_xticklabels() + S_plt.get_xticklabels()
#pylab.setp(xticklabels, visible=False)


pylab.savefig(sys.argv[5])
