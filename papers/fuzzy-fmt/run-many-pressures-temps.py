#!/usr/bin/env python2
from __future__ import division
import shlex, subprocess, os
import numpy as np
########################################################################
## This script is used to run new-soft
########################################################################

## System Size Determination
spheres = 108 # Number of Spheres
# Simulation Characteristics
iters = 1e4

# Densities to test
rho = [0.98,1.2,100.9]

# Temperatures to test
T = [1.0]

# Directory and Filename Information
directory = "data"

try:
	os.makedirs(directory)
except:
	pass
for density in rho:
	for temp in T:
		dr = 0.01*(float(spheres)/density)**(1/3)
		filename = "ff-"+str(density)+"_temp-"+str(temp)
        #~ command_line = "../../new-soft "\
		command_line = "echo rq run -J {filename} ../../new-soft "\
		"--density {density} --temp {temp} --sphereNum {spheres}"\
		" --iters {iters:.0f} --dr {dr} --dir {directory} "\
		"--filename {filename}".format(
		                  density=density, 
		                  temp=temp, 
		                  spheres=spheres,
						  iters=iters,
						  dr=dr,
						  directory=directory,
						  filename=filename,)
		print command_line
		args = shlex.split(command_line)
		p = subprocess.call(args)
