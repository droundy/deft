#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) != 7:
    print "usage:  python2 %s ww ff lenx lenyz min_T iterations" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

assert os.path.exists('SConstruct')
# switch to deft project directory and build SWMC
assert not os.system('fac square-well-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
lenx = float(sys.argv[3])
lenyz = float(sys.argv[4])
min_T = float(sys.argv[5])
iterations = round(float(sys.argv[6]))

n = ff/(4*numpy.pi/3)
N = round(n*lenyz*lenyz*lenx)
mem_estimate = 10 + 0.15*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/lv'
fname = 'ww%.2f-ff%.2f-%gx%g' % (ww, ff, lenx, lenyz)

os.system('mkdir -p ' + datadir)

if os.system('which srun'):
    # srun doesn't exist so just run on this machine
    cmd = "time nice -19 ./square-well-monte-carlo"
else:
    cmd = "srun --mem=%d -J lv:%s time nice -19 ./square-well-monte-carlo" % (mem_estimate, fname)

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += ' --lenz %g --leny %g --lenx %g --sticky-wall --walls 1' % (lenyz, lenyz, lenx)

cmd += " --iterations %d --init_iters %d --golden" % (iterations, iterations)

cmd += ' --min_T %g --translation_scale 0.05' % min_T

cmd += ' --data_dir %s --filename %s' % (datadir, fname)

cmd += " > %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
os.system(cmd)
