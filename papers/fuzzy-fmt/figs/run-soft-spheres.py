#!/usr/bin/python2

from __future__ import division
from os import system
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'

srun = lambda ff, kT: 'srun -p lomem --mem=600 -J ss-%.4f-%.4f time nice -19 ' % (ff, kT)

if os.getcwd()[:5] != '/home' or os.system('srun true') != 0:
    # We are definitely not running on the cluster
    srun = lambda ff, kT: 'time nice -19 '

# always remember to build the executable before running it
system('scons -U')

def run_soft_sphere(fillingfraction, temperature):
    nspheres = round(fillingfraction*30**3/(4*pi/3))
    outfilename = 'ss-%.4f-%.4f.out' % (fillingfraction, temperature)
    system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
           (srun(fillingfraction, temperature), figsdir, fillingfraction, temperature, outfilename))

for kT in [0.1, 0.01]:
    for ff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]:
        run_soft_sphere(ff, kT)
