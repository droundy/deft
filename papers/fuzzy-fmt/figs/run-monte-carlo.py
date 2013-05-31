from __future__ import division
from os import system
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs/'
    bindir = '../..'

def run_homogeneous(fillingfraction, temperature):
    nspheres = round(fillingfraction*30**3/(4*pi/3))
    filename = '%s/mc-%.4f-%.4f' % (figsdir, fillingfraction, temperature)
    system("srun --mem=60 -J soft-%.4f-%.4f time %s/soft-monte-carlo %d 0 0.001 %s.dat periodxyz 30 kT %g > %s.out 2>&1 &" %
           (fillingfraction, temperature, bindir, nspheres, filename, temperature, filename))

def run_walls(fillingfraction, nspheres, temperature):
    if nspheres == 0:
        nspheres = round(fillingfraction*30**2*32/(4*pi/3)) # reasonable guess
    filename = '%s/mcwalls-%.4f-%.4f-%d' % (figsdir, fillingfraction, temperature, nspheres)
    system("srun --mem=60 -J softwalls-%.4f-%.4f time %s/soft-monte-carlo %d 0.01 0.001 %s.dat periodxy 30 wallz 30 kT %g > %s.out 2>&1 &" %
           (fillingfraction, temperature, bindir, nspheres, filename, temperature, filename))

run_walls(0.1, 0, 0.1)
run_walls(0.2, 0, 0.1)
run_walls(0.3, 0, 0.1)
run_walls(0.4, 0, 0.1)
run_walls(0.5, 0, 0.1)
run_walls(0.6, 0, 0.1)
run_walls(0.7, 0, 0.1)
run_walls(0.8, 0, 0.1)

run_homogeneous(0.1, 0.0001)
run_homogeneous(0.2, 0.0001)
run_homogeneous(0.3, 0.0001)
run_homogeneous(0.4, 0.0001)
run_homogeneous(0.5, 0.0001)
run_homogeneous(0.6, 0.0001)
run_homogeneous(0.7, 0.0001)
run_homogeneous(0.8, 0.0001)

run_homogeneous(0.1, 0.001)
run_homogeneous(0.2, 0.001)
run_homogeneous(0.3, 0.001)
run_homogeneous(0.4, 0.001)
run_homogeneous(0.5, 0.001)
run_homogeneous(0.6, 0.001)
run_homogeneous(0.7, 0.001)
run_homogeneous(0.8, 0.001)

run_homogeneous(0.01, 0.01)
run_homogeneous(0.1, 0.01)
run_homogeneous(0.2, 0.01)
run_homogeneous(0.3, 0.01)
run_homogeneous(0.4, 0.01)
run_homogeneous(0.5, 0.01)
run_homogeneous(0.6, 0.01)
run_homogeneous(0.7, 0.01)
run_homogeneous(0.8, 0.01)

run_homogeneous(0.1, 0.1)
run_homogeneous(0.2, 0.1)
run_homogeneous(0.3, 0.1)
run_homogeneous(0.4, 0.1)
run_homogeneous(0.5, 0.1)
run_homogeneous(0.6, 0.1)
run_homogeneous(0.7, 0.1)
run_homogeneous(0.8, 0.1)
