#!/usr/bin/python2

from __future__ import division
from os import system
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'

srun = lambda n_reduced, kT: 'srun --mem=10000 -J wca-%.4f-%.4f time nice -19 ' % (n_reduced, kT)

if os.getcwd()[:5] != '/home' or os.system('srun true') != 0:
    # We are definitely not running on the cluster
    srun = lambda n_reduced, kT: 'time nice -19 '

# always remember to build the executable before running it
system('scons -U')

def run_soft_sphere(reduced_density, temperature):
    nspheres = round(reduced_density*2**(-5.0/2.0)*30**3)
    outfilename = 'wca-%.4f-%.4f.out' % (reduced_density, temperature)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    system("%s %s/radial-wca.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature), figsdir, reduced_density, temperature, outfilename))

for kT in [0.1,.5,1,1.5,10]:
    for n_reduced in [0.2]:
        run_soft_sphere(n_reduced, kT)
