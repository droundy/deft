#!/usr/bin/env python2

#SBATCH --nodes=1 # Number of cores 
#SBATCH --ntasks=1 # Ensure that all cores are on one machine. 
#SBATCH --time=0-00:00 # Runtime in D-HH:MM. Zero is no imposed time limit.
#SBATCH --mem=100 # Memory pool for all cores (see also --mem-per-cpu).
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL. 
#SBATCH --mail-user=pommerjo@oregonstate.edu # Email to which notifications will be sent.   

from __future__ import division
import os, numpy, sys, re

if len(sys.argv) < 6:
    print "usage:  python2 %s ww ff lenx lenyz min_T method" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

if os.path.exists('run-sticky-wall.py'):
    os.chdir('../../..')

assert os.path.exists('papers')
# switch to deft project directory and build SWMC
assert not os.system('fac liquid-vapor-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
lenx = float(sys.argv[3])
lenyz = float(sys.argv[4])
min_T = float(sys.argv[5])
method_name = sys.argv[6]

if method_name == 'wltmmc':
    wl_factor = 1.0
    wl_threshold = float(sys.argv[7])
    wl_cutoff = float(sys.argv[8])
    min_important_energy = float(sys.argv[9])
    max_entropy_energy = float(sys.argv[10])

if method_name == 'samc':
    sa_t0 = float(sys.argv[7])

seed = 0
try:
    if method_name == 'wltmmc':
        seed = int(sys.argv[11])
    elif method_name == 'samc':
        seed = int(sys.argv[8])
    else:
        seed = int(sys.argv[7])
except:
    print('defaulting to seed=0')

method = ' --' + method_name
if method[-1] in "123":
    method = method[:-1] + ' --tmi-version=' + method[-1] + ' --sad-fraction=' + method[-1]

n = ff/(4*numpy.pi/3)
N = round(n*lenyz*lenyz*lenx)
mem_estimate = 10 + 0.15*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/histogram/data/lv'
if 'golden' in sys.argv:
    suffix = 'golden'
else:
    suffix = method_name
fname = 'ww%.2f-ff%.2f-%gx%g-%s' % (ww, ff, lenx, lenyz, suffix)
if method_name == 'wltmmc':
    fname += '-%g-%g' % (wl_threshold, wl_cutoff)
elif method_name == 'samc':
    fname += '-%g' % (sa_t0)
if seed != 0:
    fname += '-s%d' % seed

os.system('mkdir -p ' + datadir)

with open('%s/%s.out' % (datadir, fname), 'w') as f:
    f.write('# python %s\n' % (' '.join(sys.argv)))

cmd = "rq run -o %s/%s.out -J lv/%s ./liquid-vapor-monte-carlo" % (datadir, fname, fname)

cmd += ' --movies' # generate movie data

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

# added by JP for wltmmc 2017.
if method_name == 'wltmmc':
  cmd += " --wl-factor %g --wl-threshold %g --wl-cutoff %g" % (wl_factor, wl_threshold, wl_cutoff)
  cmd += " --min-important-energy %g --max-entropy-energy %g" % (min_important_energy, max_entropy_energy)
elif method_name == 'samc':
  cmd += ' --sa-t0 %g' % sa_t0

cmd += ' --lenz %g --leny %g --lenx %g --sticky-wall --walls 1' % (lenyz, lenyz, lenx)

cmd += " --min-samples 100000"

cmd += method

cmd += ' --min-T %g --translation-scale 0.05' % min_T

cmd += ' --dir %s --filename %s' % (datadir, fname)

cmd += ' --seed=%d' % seed

# cmd += " >> %s/%s.out 2>&1 &" % (datadir, fname)

print(cmd)
# os.system('echo %s > %s/%s.out' % (cmd, datadir, fname))
os.system(cmd)
