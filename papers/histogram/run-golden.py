from __future__ import division
import os, numpy, sys

if len(sys.argv) != 3:
    print "usage:  python %s N ff" % sys.argv[0]
    exit(1)

paperdir = 'papers/histogram'

# always remember to build the executable before running it
os.system('fac square-well-monte-carlo')

def run_golden(N, eta, ww, method):
    filename = 'golden-%s-%d-%.4f-%.4f' % (method, N, eta, ww)
    iterations = 10000000000
    cmd = ("srun --mem=600 -J %s time nice -19 ./square-well-monte-carlo --N %d --filename_suffix golden --min_T 0.1 --min_samples 10000 --%s --ww %g --iterations %d  > %s.out 2>&1 &" %
           (filename, N, method, ww, iterations, paperdir+'/data/'+filename))
    print(cmd)
    os.system(cmd)

num = int(sys.argv[1])
ff = float(sys.argv[2])

run_golden(num, ff, 1.3, 'tmmc')
