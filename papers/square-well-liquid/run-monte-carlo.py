#!/usr/bin/python2
from __future__ import division
from numpy import *
import os, sys, argparse, socket
import subprocess as sp

import arguments

# parse arguments
parser = argparse.ArgumentParser(
  description='Run monte carlo simulation(s) using sbatch.',
  parents = [arguments.parser])

args = parser.parse_args()

# useful directories
thisdir = os.path.dirname(os.path.realpath(__file__))
jobdir = thisdir+'/jobs'

# Assumes this script is placed in [deft]/papers/square-well-liquid/
projectdir = os.path.realpath(thisdir+'../../..')

# build monte carlo code
simname = 'square-well-monte-carlo'
exitStatus = sp.call(["scons","-C",projectdir,simname])
if exitStatus != 0:
    print "Build failed"
    exit(exitStatus)

for ff in args.ff:
    for interaction_scale in args.interaction_scale:
        memory = args.N # fixme: better guess
        jobname = "sw-%iN-%iw-%4.2fff" %(args.N,args.walls,ff)
        basename = "%s/%s" %(jobdir, jobname)
        scriptname = basename + '.sh'
        outname = basename + '.out'
        errname = basename + '.err'

        command = "time nice -19 %s/%s" %(projectdir, simname)

        script = open(scriptname,'w')
        script.write("#!/bin/bash\n")
        script.write("#SBATCH --mem-per-cpu=%i\n" % memory)
        script.write("#SBATCH --output %s\n" % outname)
        script.write("#SBATCH --error %s\n\n" % errname)
        script.write("echo \"Starting job with ID: %s, "
                    "Estimated memory use: %i MB.\"\n\n" %(jobname,memory))
        script.write("cd %s\n" %projectdir)
        script.write(command)
        for (arg,val) in [ ("N",args.N), ("ff",ff),("walls",args.walls),
                           ("interaction_scale",interaction_scale),
                           ("initialize",args.initialize),
                           ("iterations",args.iterations)]:
            script.write(" \\\n --" + arg + "=" + str(val))
        script.close()


        if socket.gethostname() == 'MAPHost': # if on my own computer
            sp.Popen(["bash", scriptname])
        else:
            sp.Popen(["sbatch", "-J", jobname, scriptname])
