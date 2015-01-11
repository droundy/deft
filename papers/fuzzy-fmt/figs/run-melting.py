#!/usr/bin/python2

from __future__ import division
from os import system
import os
import numpy as np
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'

srun = lambda a_reduced, n_reduced, kT: 'srun -J melting-%.4f-%.4f-%.4f time nice -19 ' % (a_reduced, n_reduced, kT)

if os.getcwd()[:5] != '/home' or os.system('srun true') != 0:
    # We are definitely not running on the cluster
    srun = lambda a_reduced, n_reduced, kT: 'time nice -19 '

# always remember to build the executable before running it
system('scons -U')

def run_melting(reduced_lattice_const, reduced_density, temperature):
    outfilename = 'melting-%.4f-%.4f-%.4f.out' % (reduced_lattice_const, reduced_density, temperature)
    #system("%s %s/new-melting.mkdat %g %g %g > %s 2>&1 &" %
    system("%s %s/new-melting.mkdat %g %g %g > %s >&1 &" %
           (srun(reduced_lattice_const, reduced_density, temperature), figsdir, reduced_lattice_const, reduced_density, temperature, outfilename))

for a_reduced in np.arange(0.8, 2.01, 0.1):
    run_melting(a_reduced, 1.4, 1.0)

