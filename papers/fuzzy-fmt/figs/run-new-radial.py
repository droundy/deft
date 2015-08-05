#!/usr/bin/python2

from __future__ import division
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'

def srun(n_reduced, kT):
    return 'srun --mem=10000 -J name-%06.4f-%04.2f nice -19' % (kT, n_reduced)

# always remember to build the executable before running it
assert not os.system('scons papers/fuzzy-fmt/figs/new-radial-lj.mkdat papers/fuzzy-fmt/figs/new-radial-wca.mkdat')

def runme(reduced_density, temperature):
    outfilename = figsdir+'/new-data/radial-lj-%06.4f-%04.2f.out' % (temperature, reduced_density)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    os.system("%s %s/new-radial-lj.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature), figsdir, reduced_density, temperature, outfilename))

def run_wca(reduced_density, temperature):
    outfilename = figsdir+'/new-data/radial-wca-%06.4f-%04.2f.out' % (temperature, reduced_density)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    cmd = ("%s %s/new-radial-wca.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature).replace('name','new-wca'), figsdir, reduced_density, temperature, outfilename))
    print(cmd)
    os.system(cmd)

def run_bh_wca(reduced_density, temperature):
    outfilename = figsdir+'/new-data/radial-bh-wca-%06.4f-%04.2f.out' % (temperature, reduced_density)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    cmd = ("%s %s/new-bh-radial-wca.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature).replace('name','new-bh-wca'), figsdir, reduced_density, temperature, outfilename))
    print(cmd)
    os.system(cmd)

def run_bh_lj(reduced_density, temperature):
    outfilename = figsdir+'/new-data/radial-bh-lj-%06.4f-%04.2f.out' % (temperature, reduced_density)
    #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
    cmd = ("%s %s/new-bh-radial-lj.mkdat %g %g > %s 2>&1 &" %
           (srun(reduced_density, temperature).replace('name','new-bh-lj'), figsdir, reduced_density, temperature, outfilename))
    print(cmd)

# runme(0.83890, 0.71)

# runme(0.957, 2.48)

# runme(0.5844, 1.235)

# runme(1.095, 2.48)
# run_wca(0.6, 10.0)
# run_wca(0.6, 5.0)
# run_wca(0.6, 2.5)

#for temp in [2.5]:
#    run_bh_wca(0.6, temp)
#    run_bh_wca(1.0, temp)
    
#run_bh_wca(0.6, 2.5)
#run_bh_wca(0.6, 0.1)

run_bh_lj(0.83890, 0.71)
run_bh_lj(0.957, 2.48)
run_bh_lj(0.5844, 1.235)

