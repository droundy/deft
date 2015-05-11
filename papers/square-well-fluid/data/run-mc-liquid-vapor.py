#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) != 6:
    print "usage:  python2 %s ww ff lenxy lenz" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

assert os.path.exists('SConstruct')
# switch to deft project directory and build SWMC
assert not os.system('fac square-well-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
lenxy = float(sys.argv[3])
lenz = float(sys.argv[4])
iterations = round(float(sys.argv[5]))

min_T = 0.1
n = ff/(4*numpy.pi/3)
N = round(n*lenxy*lenxy*lenz)
mem_estimate = 10 + 0.1*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/lv'
fname = 'ww%.2f-ff%.2f-%gx%g' % (ww, ff, lenxy, lenz)

os.system('mkdir -p ' + datadir)

if os.system('which srun'):
    # srun doesn't exist so just run on this machine
    cmd = "time nice -19 ./square-well-monte-carlo"
else:
    cmd = "srun --mem=%d -J %s time nice -19 ./square-well-monte-carlo" % (mem_estimate, fname)

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += ' --sticky-wall --lenx %g --leny %g --lenz %g' % (lenxy, lenxy, lenz)

cmd += " --iterations %d --init_iters %d --golden" % (iterations, iterations*10)

cmd += ' --de_g 0.01' # nice high-resolution radial distribution function data

cmd += ' --min_T %g' % min_T

cmd += ' --data_dir %s --filename %s' % (datadir, fname)

cmd += " > %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
os.system(cmd)
