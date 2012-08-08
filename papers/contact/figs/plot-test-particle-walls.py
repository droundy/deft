
#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, math

if len(sys.argv) != 3:
    print("Usage:  " + sys.argv[0] + " wb-filename.dat out-filename.pdf")
    exit(1)

dftdata = numpy.loadtxt(sys.argv[1])

dft_len = len(dftdata[:,0])
axis_len = math.sqrt(dft_len)
dft_dr = dftdata[2,0] - dftdata[1,0]
n0 = dftdata[:,6]
nA = dftdata[:,8]

#pylab.figure(figsize=(8, 8))
#pylab.subplots_adjust(hspace=0.001)

n_plt = pylab.subplot(1,1,1)
n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"r+-",label='wallsDFT $n$')
n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")
pylab.xlim(0,12)
pylab.ylim(-0.05,0.2)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
#     pylab.ylim(0.1,0.8)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
#     pylab.ylim(0.17,0.43)
# if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
#     pylab.ylim(0.04,0.19)
# pylab.ylim(2.9,3.1)
pylab.xlabel("position")

pylab.ylabel("filling fraction")

pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

pylab.savefig(sys.argv[2])

pylab.show()
