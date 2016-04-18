#!/usr/bin/python

from __future__ import division
from math import pi       # REALLY don't need all of math
import os, numpy as np, sys
import time
import shutil

#assert(not os.system("fac ../../../free-energy-monte-carlo ../../../free-energy-monte-carlo-infinite-case"))


R=1

i = 0
#RG recursion level

ww = 1.5

L = float((sys.argv[1]))
#arg L = [5.0]

maxN=int((L**3.0)*0.5/(4/3.0*np.pi))+1
#maxN=2
timeData=[]

N=maxN
dirname = 'scrunched-ww%4.2f-L%04.2f/i%01d/N%03d/absolute' % (ww,L,i,N)
os.system('mkdir -p '+dirname)
ff=4/3.0*np.pi*N/L/L/L+1e-9
#cmd = ' ../../../free-energy-monte-carlo --'
cmd = '../../../free-energy-monte-carlo --pack_and_save_cube'
cmd += ' --ff %g' % ff
cmd += ' --N %d' % N
cmd += ' --data_dir %s' % dirname
cmd += ' > %s/%s.out' % (dirname, "readThis")
start=time.time()
for dummy in range(0,20):#make sure first set of atoms is mixed
	os.system(cmd)
end=time.time()
timeData.append(end-start)
print "maxN=",maxN," n=",N," dt=",end-start
prevAtoms=dirname+"/baseAtoms.txt"

for N in range(maxN-1,2,-1):
    dirname = 'scrunched-ww%4.2f-L%04.2f/i%01d/N%03d/absolute' % (ww,L,i,N)
    os.system('mkdir -p '+dirname)
    shutil.copy(prevAtoms,dirname+"/baseAtoms.txt")
    prevAtoms=dirname+"/baseAtoms.txt"
    ff=4/3.0*np.pi*N/L/L/L+1e-9
    #cmd = ' ../../../free-energy-monte-carlo --'
    cmd = '../../../free-energy-monte-carlo --pack_and_save_cube'
    cmd += ' --ff %g' % ff
    cmd += ' --N %d' % N
    cmd += ' --data_dir %s' % dirname
    cmd += ' > %s/%s.out' % (dirname, "readThis")
    start=time.time()
    os.system(cmd)
    end=time.time()
    timeData.append(end-start)
    print "maxN=",maxN," n=",N," dt=",end-start
