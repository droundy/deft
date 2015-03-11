#!/usr/bin/python2

from __future__ import division
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'

def srun(n_reduced, kT):
    return 'srun --mem=10000 -J lj-%06.4f-%04.2f' % (kT, n_reduced)

# always remember to build the executable before running it
os.system('fac papers/fuzzy-fmt/figs/new-radial-lj.mkdat')

def runme(reduced_density, temperature):
    nspheres = round(reduced_density*2**(-5.0/2.0)*30**3)
    outfilename = figsdir+'/new-data/radial-lj-%06.4f-%04.2f.out' % (temperature, reduced_density)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    os.system("%s %s/new-radial-lj.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature), figsdir, reduced_density, temperature, outfilename))

runme(0.83890, 0.71)

runme(0.957, 2.48)

runme(0.5844, 1.235)

runme(1.095, 2.48)
