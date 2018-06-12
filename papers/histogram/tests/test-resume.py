#!/usr/bin/env python2

from __future__ import division, print_function
import os, numpy, sys, re

# RUN FROM DEFT
# python papers/histogram/data/run-ising-model.py 10 10E3 40E3 sad

if len(sys.argv) < 4:
    print("usage:  python3 %s N test-moves total-moves method" % sys.argv[0])
    exit(1)

#if os.path.exists('paper.tex'):
#    os.chdir('../..')

assert os.path.exists('../data')

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

datadir = 'results'
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

filenamebase = '%s-%d-%g-%g' % (method_name, N, test_moves, total_moves)

isingdir = '../figs'
cmd = ("%s/ising.exe --dir %s --N=%i --total-moves=%g --filename=%s-complete %s" 
       % (isingdir, datadir, N, total_moves, filenamebase, method))
cmd += (" && %s/ising.exe --dir %s --N=%i --total-moves=%g --filename=%s-test %s" 
        % (isingdir, datadir, N, test_moves, filenamebase, method))
cmd += (" && %s/ising.exe --dir %s --N=%i --total-moves=%g --filename=%s-test %s --resume" 
        % (isingdir, datadir, N, total_moves - test_moves, filenamebase, method))

#if method_name == 'samc':
#  cmd += ' --sa-t0 %g' % sa_t0


print(cmd)
# os.system('echo %s > %s/%s.out' % (cmd, datadir, fname))
assert(os.system(cmd) == 0)
print("checking that both results are identical...")
assert(os.system('diff -u %s/%s-complete.dat %s/%s-test.dat'
                 % (datadir, filenamebase, datadir, filenamebase)) == 0)

print("\nAll is good!!!")
