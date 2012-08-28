
#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys, math

if len(sys.argv) != 4:
    print("Usage:  " + sys.argv[0] + " wb-test-filename.dat wb-filename.dat out-filename.pdf")
    exit(1)

testdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])


dft_len = len(testdata[:,0])
axis_len = math.sqrt(dft_len)
dft_dr = testdata[2,0] - testdata[1,0]
#n0 = testdata[:,6]
#nA = testdata[:,8]
#nAdft = dftdata[:,8]

dft_len = len(dftdata[:,0])
dft_dr = dftdata[2,0] - dftdata[1,0]

stop_here = int(dft_len - 1/dft_dr)
print stop_here
start_here = int(2.5/dft_dr)
off = 1
me = 30


n_plt = pylab.subplot(1,1,1)
n_plt.plot(testdata[:,0],testdata[:,1]*4*numpy.pi/3,"r+-",label='test-particle-DFT $n$')
#n_plt.plot(testdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
#n_plt.plot(testdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")

print "dft values:"
print dftdata[2,0]
print dftdata[2,8]
print dftdata[2,5]

n_plt.plot(dftdata[:,0],dftdata[:,8]*dftdata[:,5]*4*numpy.pi/3,"g+-",markevery=me,label="$n-contact_\sigma^A$ (White Bear)")

pylab.xlim(0,12)
#pylab.ylim(-.05,0.2)
pylab.xlabel("position")

pylab.ylabel("density")

pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

pylab.savefig(sys.argv[3])

pylab.show()
