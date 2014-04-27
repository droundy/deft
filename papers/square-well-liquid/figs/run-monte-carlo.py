#!/usr/bin/python2

import os, sys, socket
import subprocess as sp

if len(sys.argv) != 5:
    print 'useage: %s ww ff N version' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])

ff = float(sys.argv[2])

N = int(sys.argv[3])

version = str(sys.argv[4])

cores = 2 if socket.gethostname() == 'MAPHost' else 8

initialize = 10000
iterations = 10000

# define some directories
figdir = os.path.dirname(os.path.realpath(__file__))
swdir = os.path.dirname(figdir)
projectdir = os.path.realpath(swdir+'../../..')
jobdir = swdir+'/jobs'
datadir = swdir+'/data/'
simname = 'square-well-monte-carlo'

# build monte carlo code
exitStatus = sp.call(["scons","-j%i"%cores,"-C",projectdir,simname])
if exitStatus != 0:
    print "Build failed"
    exit(exitStatus)

# write script to run simulation
memory = N # fixme: better guess
jobname = 'periodic-ww%02.0f-ff%02.0f-N%i-%s' \
  % (ww*100, ff*100, N, version)
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
for (arg,val) in [ ("ww",ww), ("ff",ff), ("N",N),
                   ("initialize",initialize),
                   ("iterations",iterations) ]:
    script.write(" \\\n --%s %s" %(arg,str(val)))
script.write(" \\\n --%s"%version)
script.close()

# start simulation
if socket.gethostname() == 'MAPHost':
#    sp.Popen(["bash", scriptname])
    None
else:
    sp.Popen(["sbatch", "-J", jobname, scriptname])
