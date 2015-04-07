from __future__ import division
import os, numpy

paperdir = 'papers/histogram'

# always remember to build the executable before running it
os.system('fac square-well-monte-carlo')

def run_liquid_vapor(eta, ww):
    n = eta/(4*numpy.pi/3)
    lenx = 10
    leny = 10
    lenz = 50
    nspheres = round(n*lenx*leny*lenz)
    lenx = numpy.sqrt(nspheres/n/lenz) # make mean density exact
    leny = lenx
    #width = (nspheres/density)**(1.0/3) # to get density just right!
    filename = 'lv-%.4f-%.4f' % (eta, ww)
    cmd = ("srun --mem=600 -J %s time nice -19 ./square-well-monte-carlo --filename %s --N %d --lenx %g --leny %g --lenz %g --tmmc --ww %g --ff 0 --sticky-wall > %s.out 2>&1 &" %
           (filename, filename, nspheres, lenx, leny, lenz, ww, paperdir+'/data/'+filename))
    print(cmd)
    os.system(cmd)

run_liquid_vapor(0.2, 1.5)
