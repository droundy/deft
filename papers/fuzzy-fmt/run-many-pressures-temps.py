#!/usr/bin/env python2

import shlex, subprocess, os
import numpy as np

# Thermodynamic characteristics
fillFrac = 0.8
temp = 0.2

## System Size Determination
lenx = 10
leny = 10
lenz = 10
wallx = 0
wally = 0
wallz = 0

# Simulation Characteristics
iters = 1000000
dr = 0.1

# Directory and Filename Information
direc_name = "CHRIS"


ffmin = 0.1
ffmax = 1.1
dff = 0.1

tempmin = 0.01
tempmax = 3.5
dtemp = 0.5


tempFF = np.arange(ffmin,ffmax,dff)
tempTemp = np.arange(tempmin,tempmax,dtemp)

try:
	os.makedirs('data')
except:
	pass

for fillFrac in tempFF:
	for temp in tempTemp:
		filename = "ff-"+str(fillFrac)+"_temp-"+str(temp)
		command_line = "sbatch -o data/{}.out -c 1 -n 1 -J {} ./run-new-soft.py --density {} --temp {}".format(
		                  filename, filename, fillFrac, temp)
		args = shlex.split(command_line)
		
		#print args
		p = subprocess.call(args)




