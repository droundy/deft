#!/usr/bin/env python2
import shlex, subprocess, os
import numpy as np
########################################################################
## This script is used to run a many instances of new-soft
########################################################################

# What values do we tell the system to use.
NDT = 1		# Number, Density, Temperature. System calcs Volume.
NVT = 0		# Number, Volume, Temperature.	You calc Density.
DVT = 0		# Density, Volume, Temperature.	System calcs Number
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
# Simulation Characteristics
iters = 10000000
dr = 0.1
# Densities to test
pmin = 0.1
pmax = 1.1
dp = 0.1
# Temperatures to test
Tmin = 0.01
Tmax = 3.5
dT = 0.5
# Makes a new directory for you if it's not already there.
try:
	os.makedirs('data')
except:
	pass

for density in np.arange(pmin,pmax,dp):
	for temp in np.arange(Tmin,Tmax,dT):
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
