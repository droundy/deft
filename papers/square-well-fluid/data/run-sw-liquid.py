#!/usr/bin/python2

from __future__ import division
from os import system
import os
from math import pi

figsdir = 'papers/square-well-fluid/figs/'

def srun(density, lam, kT):
    return 'srun --mem=5000 -J sw-fluid-%.2f-%.2f-%.2f time nice -19 ' % (density, lam,  kT)

if os.getcwd()[:5] != '/home' or os.system('srun true') != 0:
    # We are definitely not running on the cluster
    srun = lambda density, lam,  kT: 'time nice -19 '

# always remember to build the executable before running it
assert not system('scons papers/square-well-fluid/figs/new-radial-sw.mkdat')

def run_sw_fluid(density, lam, temperature):
    outfilename = 'sw-fluid-%.2f-%.2f-%.2f.out' % (density, lam, temperature)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    system("%s %s/new-radial-sw.mkdat %g %g %g > %s 2>&1 &" %
           (srun(density, lam, temperature), figsdir, density, lam, temperature, outfilename))

for kT in [1, 2.0, 3.0, 4.0, 5.0, 10.0]:
    for density in [0.2]:
        run_sw_fluid(density, 1.3, kT)
