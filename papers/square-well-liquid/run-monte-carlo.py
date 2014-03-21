#!/usr/bin/python2
from numpy import *
import os, sys, argparse, socket
import subprocess as sp

import arguments
from paramID import *

# parse arguments
parser = argparse.ArgumentParser(
  description='Run monte carlo simulation(s) using sbatch.',
  parents = [arguments.parser])

parser.add_argument(
		'-initialize', metavar='INT', type=int, default=1000000,
		help='Number of iterations to run for initialization')

parser.add_argument(
		'-iterations', metavar='INT', type=int, default=1000000000000,
		help='Number of simulation iterations')

args = parser.parse_args()

# build monte carlo code
exitStatus = sp.call(["scons","-C",projectdir,simname])
if exitStatus != 0:
    print "Build failed"
    exit(exitStatus)

# store all sets to run
paramList = []
for ww in args.ww:
    for ff in args.ff:
        for N in args.N:
            paramList.append(
                paramID(args.walls,ww,ff,N,args.weights))

for p in paramList:
    memory = p.N # fixme: better guess
    jobname = p.name('N')
    basename = "%s/sw-%s" %(jobdir, jobname)
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
                 "Estimated memory use: %i MB.\"\n\n"
                 %(jobname,memory))
    script.write("cd %s\n" %projectdir)
    script.write(command)
    for (arg,val) in [ ("ww",p.ww), ("ff",p.ff), ("N",p.N),
                       ("walls",args.walls),
                       ("initialize",args.initialize),
                       ("iterations",args.iterations) ]:
        script.write(" \\\n --" + arg + "=" + str(val))
    if p.weights: script.write(" \\\n --weights")
    script.close()

    if socket.gethostname() == 'MAPHost':
        sp.Popen(["bash", scriptname])
    else:
        sp.Popen(["sbatch", "-J", jobname, scriptname])
