#!/usr/bin/python2

from __future__ import division
import os, numpy, sys

i = eval(sys.argv[1])
#RG recursion level

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float((sys.argv[3]))
#arg L = [5.0]

L_i = L*(2**i)
#incorperation of RG parameter i

Ns = eval(sys.argv[4])
#arg Ns = [range(2,10)]

Overwrite = False
if '-O' in sys.argv:  # Check for overwrite flag in arguments
    Overwrite = True

for N in Ns:
    dirname = 'scrunched-ww%4.2f-L%04.2f/i%01d/N%03d' % (ww,L,i,N)
    if Overwrite:
        os.system('rm -rf ' + dirname)
    print('mkdir -p ' + dirname)
    os.system('mkdir -p '+ dirname)
    filename = 'ww%4.2f-L%04.2f-N%03d' % ( ww, L, N)
    iterations = 1000000


    #cmd = 'srun -J %s' % filename   
    cmd = ' ../../../square-well-monte-carlo'
    cmd += ' --filename %s' % filename
    cmd += ' --N %d' % N
    cmd += ' --min_T 0.1'
    cmd += ' --data_dir %s' % dirname
    cmd += ' --min_samples 10000'
    cmd += ' --golden'
    cmd += ' --ww %g' % ww
    cmd += ' --iterations %d' %iterations
    cmd += ' --lenx %g --leny %g --lenz %g' % (L_i,L_i,L_i) 
    cmd += ' > %s/%s.out 2>&1 &' % (dirname, filename) 
    # cmd = ("../../../square-well-monte-carlo --filename %s --N %d --min_T 0.1 --min_samples 10000 --golden --ww %g --iterations %d --lenx %g --leny %g --lenz %g --data_dir .  > %s.out 2>&1 &" %
    #        (filename, N, ww, iterations, L, L, L, filename))
    print(cmd)
    os.system(cmd)


