#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) != 4:
    print "usage:  python2 %s translation_scale translation_max iterations" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

assert os.path.exists('SConstruct')
# switch to deft project directory and build SWMC
assert not os.system('fac square-well-monte-carlo')

ww = 1.3
ff = 0.1
lenx = 50
lenyz = 10
min_T = 0.1
translation_scale = float(sys.argv[1])
translation_max = float(sys.argv[2])
iterations = round(float(sys.argv[3]))

n = ff/(4*numpy.pi/3)
N = round(n*lenyz*lenyz*lenx)
mem_estimate = 10 + 0.1*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/timings'
fname = 'max-%g-scale-%g' % (translation_max, translation_scale)

os.system('mkdir -p ' + datadir)

if os.system('which srun'):
    # srun doesn't exist so just run on this machine
    cmd = "time ./square-well-monte-carlo"
else:
    cmd = "srun --mem=%d -J %s time ./square-well-monte-carlo" % (mem_estimate, fname)

cmd += ' --translation_scale %g --translation_max %g' % (translation_scale, translation_max)

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += ' --lenz %g --leny %g --lenx %g --sticky-wall --walls 1' % (lenyz, lenyz, lenx)

cmd += " --iterations %d --init_iters %d --tmmc" % (iterations, 10*iterations)

cmd += ' --min_T %g' % min_T

cmd += ' --data_dir %s --filename %s' % (datadir, fname)

cmd += " > %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
os.system(cmd)
