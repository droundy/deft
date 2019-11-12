# This script computes the number of iterations in the Poulaine paper

import numpy as np

lnf = 1

M = 50
N = 100

iters = 0

for main in range(0,M):
    for sub in range(0,N):
        lnf *= 0.9
        iters += 1e4/np.sqrt(lnf)
        print('   substage: ', lnf, 'gaining', 1.0/np.sqrt(lnf), '-> total', iters)
    lnf /= 0.9**(N-1)
    print('returning to', lnf)

print('finally', '%e' % iters)
