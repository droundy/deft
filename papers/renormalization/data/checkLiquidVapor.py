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
	maxPercentDone=0.0
	with open(fbase+'lv-data.out','r') as file:
		for line in file:
			if "done" in line:
				nextPercentDone=float(line.split('% done')[0].split('(')[1])
				maxPercentDone=max(maxPercentDone,nextPercentDone)
	if maxPercentDone<100.0: return False
	return True

