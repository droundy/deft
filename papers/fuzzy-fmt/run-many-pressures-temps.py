#!/usr/bin/env python2
from __future__ import division
import subprocess, os
import numpy as np
########################################################################
## This script is used to run new-soft
########################################################################

## System Size Determination
spheres = 256
 # Number of Spheres
# Simulation Characteristics
iters = 1e9

# Densities to test
rho = [1.1,1.2,0.9]

# Temperatures to test
T = [2.0]

# Directory and Filename Information
directory = "data2"

try:
	os.makedirs(directory)
except:
	pass
for density in rho:
	for temp in T:
		dr = 0.01*(float(spheres)/density)**(1/3)
		filename = "ff-"+str(density)+"_temp-"+str(temp)
        #~ command_line = "../../new-soft "\
		command_line = " rq run -J {filename} ../../new-soft "\
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
        p = subprocess.call(command_line, shell=True)
# srun -J sleep-quietly nohup nice -19 hostname &> hostname.out &
# scancel
# man srun
