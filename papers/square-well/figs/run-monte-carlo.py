#!/usr/bin/python2
from __future__ import division
from numpy import *
import os, sys, argparse
import subprocess as sp

import arguments

# parse arguments
parser = argparse.ArgumentParser(
  description='Run monte carlo simulation(s) using sbatch.',
  parents = [arguments.parser])

args = parser.parse_args()

# useful directories
thisdir = os.path.dirname(os.path.realpath(__file__))
figdir = thisdir
jobdir = figdir+'/jobs'

# Assumes this script is placed in [deft]/papers/square-well/figs/
projectdir = os.path.realpath(thisdir+'../../../..')

# build monte carlo code
simname = 'square-well-monte-carlo'
sp.call(["scons","-C",projectdir,simname])

for ff in args.ff:
    memory = args.N/40 # fixme: better guess
    jobname = "sw-%iN-%iw-%4.2fff" %(args.N,args.walls,ff)
    basename = "%s/%s" %(jobdir, jobname)
    scriptname = basename + '.tmp.sh'
    outname = basename + '.out'

    command = "time nice -19 %s/%s" %(projectdir, simname)

    script = open(scriptname,'w')
    script.write("#!/bin/bash\n")
    script.write("#SBATCH --mem-per-cpu=%i\n" % memory)
    script.write("#SBATCH --output %s\n\n" % outname)
    script.write("echo \"Starting job with ID: %s, "
                 "Estimated memory use: %i MB.\"\n\n" %(jobname,memory))
    script.write("cd %s\n" %projectdir)
    script.write(command)
    for (arg,val) in [ ("ff",ff), ("N",args.N), ("walls",args.walls),
                       ("iterations",args.iterations),
                       ("initialization_iterations",
                        args.initialization_iterations)]:
        script.write(" \\\n --" + arg + "=" + str(val))
    script.close()

    sp.Popen(["sbatch", "-J", jobname, scriptname])
