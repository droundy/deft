#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) != 4:
    print "usage:  python2 %s ww ff N" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

assert os.path.exists('SConstruct')
# switch to deft project directory and build SWMC
assert not os.system('fac square-well-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
N = int(sys.argv[3])


min_T = 0.1
iterations = 1e6
mem_estimate = 10 + 0.1*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/mc'
fname = 'ww%.2f-ff%.2f-N%d' % (ww, ff, N)

cmd = "srun --mem=%d -J %s time nice -19 ./square-well-monte-carlo" % (mem_estimate, fname)
#cmd = "time nice -19 ./square-well-monte-carlo"

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += " --iterations %d --init_iters %d --golden" % (iterations, iterations)

cmd += ' --data_dir %s --filename %s' % (datadir, fname)

cmd += " > %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
os.system(cmd)

