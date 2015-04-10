#!/usr/bin/python2

from __future__ import division
import os, numpy, sys

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[2])
#arg L = [5.0]

Ns = eval(sys.argv[3])
#arg Ns = [range(2,10)]

for N in Ns:
    filename = 'data/ww%4.2f-L%04.2f-N%03d' % (ww, L, N)
    iterations = 1000000
    cmd = ("../../square-well-monte-carlo --filename %s --N %d --min_T 0.1 --min_samples 10000 --tmmc --ww %g --iterations %d --lenx %g --leny %g --lenz %g --data_dir .  > %s.out 2>&1 &" %
           (filename, N, ww, iterations, L, L, L, filename))
    print(cmd)
    os.system(cmd)
