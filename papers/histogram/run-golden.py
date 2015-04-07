#!/usr/bin/env python2
from __future__ import division
import os, numpy, sys, re

if len(sys.argv) != 4:
    print "usage:  python %s ww ff N" % sys.argv[0]
    exit(1)

paperdir = 'papers/histogram'

# switch to deft project directory and build SWMC
filepath = os.path.abspath(__file__)
deft_dir = re.sub('deft/.*','deft',filepath)
paper_dir = re.sub('histogram/.*','histogram',filepath)
data_dir = paper_dir+'/data'
os.chdir(deft_dir)
os.system('fac square-well-monte-carlo')

def run_golden(ww, ff, N, min_T, method):
    filename = 'golden-%s-N%d-ff%.0f-ww%.0f' % (method, N, ff*100, ww*100)
    iterations = 10000000000
    min_samples = 10000
    cmd = ("srun --mem=600 -J %s time nice -19 ./square-well-monte-carlo --ww %g --ff %g --N %d --min_T %g --%s --iterations %d --min_samples %d --filename_suffix golden > %s.out 2>&1 &" %
           (filename, ww, ff, N, min_T, method, iterations, min_samples, data_dir+'/'+filename))
    print(cmd)
    os.system(cmd)

ww = float(sys.argv[1])
ff = float(sys.argv[2])
N = int(sys.argv[3])
min_T = 0.2
method = 'tmmc'

run_golden(ww, ff, N, min_T, method)
