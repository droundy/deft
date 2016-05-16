#!/usr/bin/python

from __future__ import division
from math import pi       # REALLY don't need all of math
import os, numpy as np, sys
import glob
import checkAbsolute as absolute
import listAll
def absolute_f(fbase=None):
    if fbase==None: return None
    # find the partition function yielding the absolute free energy  using 'absolute/' data
    if fbase[-1]=="/": fbase=fbase[:-1]
    fbase = fbase+ '/absolute/'
    num_files=0
    while os.path.isfile(fbase+"%05d.dat"%num_files):
		num_files+=1
    #print num_files
    #num_files = len(glob .glob(fbase+'*.dat')) # figure out how many files there are
    successes = 0
    total = 0
    ratios = np.zeros(num_files)
    absolute_f = 0
   
    #print("Num_files is: %s" % num_files)
    for j in xrange(0,num_files):
        filename = fbase + '%05d' % (j)
        #if j==num_files-1:
			#print 'filename is "%s" and j is %d' % (filename, j)
        with open(filename+".dat") as file:
            for line in file:
                #print line
                if ("N: " in line):
                    N = int(line.split()[-1])
                if("ff_small: " in line):
                    ff = float(line.split()[-1])
                #if("failed small checks: " in line):
                if("total valid small" in line):
                    successes = float(line.split()[-1])
                    #print successes
                #if("valid small checks: " in line):
                if("total checks of small cell: " in line):
                    total = float(line.split()[-1])
                    #print total
                    #break
        ratios[j] = successes/total
        absolute_f += -np.log(ratios[j])

    #print("Ratios array is: %s" % ratios)
    #print("Calculated absolute_f is: %g" % absolute_f)
    #print("Compare with: %g" % ((4*ff - 3*ff**2)/(1-ff)**2))
    return absolute_f


listAll.findFolders()
for i in range(0,len(listAll.validFoldersss)):
	for j in range(0,len(listAll.validFoldersss[i])):
		print listAll.validFoldersss[i][j][0]
		for k in range(0,len(listAll.validFoldersss[i][j])):
			if False==absolute.check(fbase=listAll.validFoldersss[i][j][k]):
				continue
			absoluteF=absolute_f(fbase=listAll.validFoldersss[i][j][k])
			print "N",listAll.validFoldersss[i][j][k].split("/N")[1].split("/")[0],absoluteF
			#print listAll.validFoldersss[i][j][k]
			with open(listAll.validFoldersss[i][j][k]+"absolute/Sexc.dat",'w') as file:
				file.write("%f"%(-absoluteF))
