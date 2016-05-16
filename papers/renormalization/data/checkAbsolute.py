#!/usr/bin/env python2
from __future__ import division
import os, sys, re
import os
import shutil
import numpy as np

#sorts the data files by:


def check(fbase=None):
	if fbase==None: return
	if len(fbase)>0 and fbase[-1]!='/': fbase=fbase+'/'
	if len(fbase)>0 and fbase[0]!='/': base='/'+fbase
	outFiles=[]
	dosFiles=[]#does files is to verify a dos file was created
	L=None
	ww=None
	if False==os.path.isfile(fbase+'lv-data.out'): return False
	with open(fbase+'lv-data.out','r') as file:
		for line in file:
			if '--L' in line: L=float(line.split('--L ')[1])
			if '--w' in line: ww=float(line.split('--ww ')[1])
	if L==None: print "checkAbsolute did not find L"; return False
	dirname = fbase+'absolute/'
	if not os.path.isdir(dirname): return False
	i=-1
	Ls=[]
	while True:
		Lnow=None
		Lscalar=None
		i+=1
		fname=dirname+'%.5d.out'%i
		if not os.path.isfile(fname): break
		read=open(fname,'r')
		lines=read.readlines()
		read.close()
		maxPercentDone=0.0
		for line in lines:
			if 'complete' in line:
				nextPercentDone=float(line.split('(')[1].split('%')[0])
				maxPercentDone=max(maxPercentDone,nextPercentDone)
			if 'cell dimensions' in line:
				Lnow=float(line.split('(')[1].split(',')[0])
			if 'scaling factor' in line:
				Lscalar=float(line.split('Using scaling factor of ')[1])
		if maxPercentDone<100.0:
			break
		if Lnow==None:
			print "could not find current box size"
			break
		if i==0:
			Ls.append([1e200,Lnow])
			continue
		if Lscalar==None:
			print "could not find scale factor"
			break
		Ls.append([Lnow,Lnow*Lscalar])
	if len(Ls)==0: return False
	notSorted=True
	while notSorted==True:
		notSorted=False#start as possibly sorted
		for i in range(0,len(Ls)-1):
			if Ls[i][0]<Ls[i+1][0]:
				notSorted=True
				Ls[i],Ls[i+1]=Ls[i+1],Ls[i]
	if abs(Ls[-1][1]-L)>0.0001: return False
	for i in range(0,len(Ls)-1):
		if abs(Ls[i][1]-Ls[i+1][0])>0.0001: return False
	return True
