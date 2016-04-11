#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) < 7:
    print "usage:  python2 %s ww ff lenx lenyz min_T iterations" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

if os.path.exists('run-mc-liquid-vapor.py'):
    os.chdir('../../..')

assert os.path.exists('papers')
# switch to deft project directory and build SWMC
assert not os.system('fac liquid-vapor-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
lenx = float(sys.argv[3])
lenyz = float(sys.argv[4])
min_T = float(sys.argv[5])

if 'tmi' in sys.argv:
    method = ' --tmi'
elif 'toe' in sys.argv:
    method = ' --toe'
else:
    method = ' --tmmc'

n = ff/(4*numpy.pi/3)
N = round(n*lenyz*lenyz*lenx)
mem_estimate = 10 + 0.15*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/lv'
fname = 'ww%.2f-ff%.2f-%gx%g-tmmc' % (ww, ff, lenx, lenyz)
if 'tmi' in sys.argv:
    fname = 'ww%.2f-ff%.2f-%gx%g-tmi' % (ww, ff, lenx, lenyz)
elif 'toe' in sys.argv:
    fname = 'ww%.2f-ff%.2f-%gx%g-toe' % (ww, ff, lenx, lenyz)

os.system('mkdir -p ' + datadir)

if os.system('which srun'):
    # srun doesn't exist so just run on this machine
    cmd = "time nice -19 ./liquid-vapor-monte-carlo"
else:
    cmd = "srun --mem=%d -J lv:%s time nice -19 ./liquid-vapor-monte-carlo" % (mem_estimate, fname)

cmd += ' --movies' # generate movie data

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += ' --lenz %g --leny %g --lenx %g --sticky-wall --walls 1' % (lenyz, lenyz, lenx)

cmd += " --min-samples 100"

cmd += method

cmd += ' --min-T %g --translation-scale 0.05' % min_T

cmd += ' --dir %s --filename %s' % (datadir, fname)

cmd += " >> %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
os.system('echo %s > %s/%s.out' % (cmd, datadir, fname))
os.system(cmd)
