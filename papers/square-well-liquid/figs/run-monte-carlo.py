#!/usr/bin/python2
import os, sys, argparse, socket
import subprocess as sp

import arguments
from info import *

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
for ww in args.ww:
    for ff in args.ff:
        for N in args.N:
            for kT in args.kT:
                memory = N # fixme: better guess
                jobname = 'periodic-ww%04.2f-ff%04.2f-N%i' % (ww, ff, N)
                if kT != 0: jobname += '-kT%g' %kT
                if args.flat: jobname += '-flat'
                elif args.nw: jobname += '-nw'
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
                             "Estimated memory use: %i MB.\"\n\n"
                             %(jobname,memory))
                script.write("cd %s\n" %projectdir)
                script.write(command)
                for (arg,val) in [ ("ww",ww), ("ff",ff), ("N",N), ("kT",kT),
                                   ("walls",args.walls),
                                   ("initialize",args.initialize),
                                   ("iterations",args.iterations) ]:
                    script.write(" \\\n --" + arg + "=" + str(val))
                if args.flat: script.write(" \\\n --flat")
                elif args.nw: script.write(" \\\n --nw")
                script.close()

                if socket.gethostname() == 'MAPHost':
                    #sp.Popen(["bash", scriptname])
                    None
                else:
                    sp.Popen(["sbatch", "-J", jobname, scriptname])
