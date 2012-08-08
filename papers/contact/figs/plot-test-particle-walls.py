
#!/usr/bin/python

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

n_plt = pylab.subplot(1,1,1)
n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"r+-",label='wallsDFT $n$')
n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")
pylab.xlim(0,12)
pylab.ylim(-0.05,0.2)
pylab.xlabel("position")

pylab.ylabel("filling fraction")

pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

pylab.savefig(sys.argv[2])

pylab.show()
