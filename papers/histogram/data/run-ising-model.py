#!/usr/bin/env python2

from __future__ import division
import os, numpy, sys, re

# RUN FROM DEFT
# python papers/histogram/data/run-ising-model.py 10 10E3 40E3 sad

if len(sys.argv) < 4:
    print "usage:  python2 %s N test-moves total-moves method" % sys.argv[0]
    exit(1)

#if os.path.exists('paper.tex'):
#    os.chdir('../..')

if os.path.exists('run-ising-model.py'):
    os.chdir('../../..')

assert os.path.exists('papers')
# switch to deft project directory and build ising.exe
assert not os.system('fac papers/histogram/figs/ising.exe')

N = int(sys.argv[1])
test_moves = float(sys.argv[2])
total_moves = float(sys.argv[3])
method_name = sys.argv[4]

if method_name == 'samc':
    sa_t0 = float(sys.argv[7])

seed = 0
try:
    if method_name == 'samc':
        seed = int(sys.argv[8])
    else:
        seed = int(sys.argv[7])
except:
    print('defaulting to seed=0')

method = ' --' + method_name

datadir = 'papers/histogram/data/ising'
if 'golden' in sys.argv:
    suffix = 'golden'
else:
    suffix = method_name
fname = 'N%i-moves%.1E-%s' % (N, total_moves, suffix) 
if method_name == 'samc':
    fname += '-%g' % (sa_t0)
if seed != 0:
    fname += '-s%d' % (seed)

os.system('mkdir -p ' + datadir) # create ising folder in directory

# create a file in datadir with fname
#with open('%s/%s.out' % (datadir, fname), 'w') as f:
#    f.write('# python %s\n' % (' '.join(sys.argv)))

isingdir = 'papers/histogram/figs'
cmd = "%s/ising.exe --N=%i --total-moves=%g --filename=ising-complete %s" % (isingdir, N, total_moves, method)
cmd += "&& %s/ising.exe --N=%i --total-moves=%g --filename=ising-test %s" % (isingdir, N, test_moves, method)
cmd += "&& %s/ising.exe --N=%i --total-moves=%g --filename=ising-test %s --resume" % (isingdir, N, total_moves - test_moves, method)

#if method_name == 'samc':
#  cmd += ' --sa-t0 %g' % sa_t0


print(cmd)
# os.system('echo %s > %s/%s.out' % (cmd, datadir, fname))
os.system(cmd)
