#!/usr/bin/env python2
import shlex, subprocess, os
import numpy as np
########################################################################
## This script is used to run new-soft
########################################################################

## System Size Determination
spheres = 108 # Number of Spheres
# Simulation Characteristics
iters = 1e9
dr = 0.1
# Densities to test
pmin = 0.1
pmax = 0.2
dp = 0.1
# Temperatures to test
Tmin = 0.01
Tmax = 0.51
dT = 0.5
# Makes a new directory for you if it's not already there.
try:
	os.makedirs('data8-10-17')
except:
	pass

# Directory and Filename Information
directory = "data8-10-17"

for density in np.arange(pmin,pmax,dp):
	for temp in np.arange(Tmin,Tmax,dT):
		filename = "ff-"+str(density)+"_temp-"+str(temp)
		#command_line = "sbatch -o data/{}.out -c 1 -n 1 -J {}"\
		command_line = "rq run -J {filename} ../../new-soft --density {} --temp {} --sphereNum {}"\
		" --iters {iters:.0f} --dr {dr} --dir {directory}".format(
		                  density, temp, spheres,
						  iters=iters,
						  dr=dr,
						  directory=directory,
						  filename=filename)
		print command_line
		args = shlex.split(command_line)
		p = subprocess.call(args)
