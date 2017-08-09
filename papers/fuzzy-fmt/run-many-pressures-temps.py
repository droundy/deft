#!/usr/bin/env python2

import shlex, subprocess, os
import numpy as np

# What values do we tell the system to use.
NDT = 1		# Number, Density, Temperature On
NVT = 0		# Number, Volume, Temperature Off
DVT = 0		# Density, Volume, Temperature Off

## System Size Determination
spheres = 108 # Number of Spheres
lenx = 10
leny = 10
lenz = 10
lenxy = 10

# Wall Existence
wallx = 0
wally = 0
wallz = 0

#

# Simulation Characteristics
iters = 10000000
dr = 0.1

ffmin = 0.1
ffmax = 1.1
dff = 0.1

tempmin = 0.01
tempmax = 3.5
dtemp = 0.5


try:
	os.makedirs('data')
except:
	pass

for density in np.arange(ffmin,ffmax,dff):
	for temp in np.arange(tempmin,tempmax,dtemp):
		filename = "ff-"+str(density)+"_temp-"+str(temp)
		command_line = "sbatch -o data/{}.out -c 1 -n 1 -J {}"\
		" ./run-new-soft.py --density {} --temp {} --number {}"\
		" --NDT {} --NVT {} --DVT {} --lenx {} --leny {} --lenz {}"\
		" --wallx {} --wally {} --wallz {}"\
		" --iters {} --dr {}".format(
		                  filename, filename,
		                  density, temp, spheres,
						  NDT, NVT, DVT, lenx, leny, lenz,
						  wallx, wally, wallz,
						  iters, dr)
		args = shlex.split(command_line)
		p = subprocess.call(args)




