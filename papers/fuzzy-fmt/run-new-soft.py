#!/usr/bin/env python2
import shlex, subprocess, sys, argparse, os
import numpy as np

########################################################################
## This script is used to run a single instance of new-soft
########################################################################

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
parser.add_argument('--NDT', metavar='SIMTYPE', type=int, action='store',
                    default=0,
                    help='Input Num,Density,Temp.')
parser.add_argument('--DVT', metavar='SIMTYPE', type=int, action='store',
                    default=0,
                    help='Input Density,Vol,Temp.')
parser.add_argument('--NVT', metavar='SIMTYPE', type=int, action='store',
                    default=0,
                    help='Input Num,Vol,Temp.')
parser.add_argument('--lenx', metavar='LENGTH', type=float, action='store',
                    default=0,
                    help='System Length in X')
parser.add_argument('--leny', metavar='LENGTH', type=float, action='store',
                    default=0,
                    help='System Length in Y.')
parser.add_argument('--lenz', metavar='LENGTH', type=float, action='store',
                    default=0,
                    help='System Length in Z')
parser.add_argument('--lenxy', metavar='LENGTH', type=float, action='store',
                    default=0,
                    help='Input Num,Vol,Temp.')
parser.add_argument('--wallx', metavar='WALL', type=int, action='store',
                    default=0,
                    help='Wall on X')
parser.add_argument('--wally', metavar='WALL', type=int, action='store',
                    default=0,
                    help='Wall on Y')
parser.add_argument('--wallz', metavar='WALL', type=int, action='store',
                    default=0,
                    help='Wall on Z')
parser.add_argument('--iters', metavar='ITERATIONS', type=int, action='store',
                    default=0,
                    help='Total Iterations')
parser.add_argument('--dr', metavar='LENGTH', type=float, action='store',
                    default=0,
                    help='Random Move Length Max')
args = parser.parse_args()

# Directory and Filename Information
direc_name = "data8-10-17"
## Make this Directory name a parsed argument

filename = "ff-"+str(args.density)+"_temp-"+str(args.temp)
command_line ="../../new-soft --NDT "+str(args.NDT)+" --DVT "+str(args.DVT)+\
				" --NVT "+str(args.NVT)+" --lenx "+str(args.lenx)+" --leny "+\
				str(args.leny)+" --lenz "+str(args.lenz)+" --lenxy "+str(args.lenxy)+\
				" --sphereNum " + str(args.number)+" --density "+\
				str(args.density)+" --temp "+str(args.temp)+" --wallx "\
				+str(args.wallx)+" --wally "+str(args.wally)+" --wallz "+\
				str(args.wallz)+" --iters "+str(args.iters)+" --dr "+str(args.dr)+\
				" --dir "+direc_name+" --filename "+filename

print command_line
p = subprocess.call(shlex.split(command_line))




