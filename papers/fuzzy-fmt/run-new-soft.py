#!/usr/bin/env python2

import shlex, subprocess, sys, argparse, os
import numpy as np

parser = argparse.ArgumentParser(description='run new-soft')

parser.add_argument('--density', metavar='DENSITY', type=float, action='store',
                    default=0.2,
                    help='the reduced density')
parser.add_argument('--temp', metavar='TEMPERATURE', type=float, action='store',
                    default=1.0,
                    help='the reduced temperature')
parser.add_argument('--number', metavar='SPHERES', type=int, action='store',
                    default=1,
                    help='the number of spheres')
args = parser.parse_args()

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
direc_name = "data"

filename = "ff-"+str(args.density)+"_temp-"+str(args.temp)
command_line ="../../new-soft --lenx "+str(lenx)+" --leny "+str(leny)+\
				" --lenz "+str(lenz)+" --fillFrac "+str(args.density)+\
				" --temp "+str(args.temp)+" --wallx "+str(wallx)+\
				" --wally "+str(wally)+" --wallz "+str(wallz)+\
				" --iters "+str(iters)+" --dr "+str(dr)+\
				" --dir "+direc_name+" --filename "+filename

print command_line
p = subprocess.call(shlex.split(command_line))




