#!/usr/bin/python2

import os, sys, socket
import subprocess as sp

if not len(sys.argv) in [5,6]:
    print 'useage: %s ww ff N versions filename_suffix' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])

ff = float(sys.argv[2])

N = int(sys.argv[3])

versions = eval(sys.argv[4])

if len(sys.argv) == 6:
    filename_suffix = sys.argv[5]

initialize = 10000
iterations = 10000

# define some directories
swdir = os.path.dirname(os.path.realpath(__file__))
figdir = os.path.realpath(swdir+'/figs')
projectdir = os.path.realpath(swdir+'/../..')
jobdir = swdir+'/jobs'
datadir = swdir+'/data'
simname = 'square-well-monte-carlo'

cores = 2 if socket.gethostname() == 'MAPHost' else 8

# build monte carlo code
exitStatus = sp.call(["scons","-j%i"%cores,"-C",projectdir,simname],
                     stdout = open(os.devnull,"w"),
                     stderr = open(os.devnull,"w"))
if exitStatus != 0:
    print "scons failed"
    exit(exitStatus)

for version in versions:
    memory = 20*N # fixme: better guess
    jobname = 'periodic-ww%04.2f-ff%04.2f-N%i-%s' %(ww, ff, N, version)
    basename = "%s/%s" %(jobdir, jobname)
    scriptname = basename + '.sh'
    outname = basename + '.out'
    errname = basename + '.err'

    command = "time %s/%s" %(projectdir, simname)

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
    script.write(" \\\n --%s"%version.replace("kT","kT "))
    try:
        script.write(" \\\n --filename_suffix %s" %filename_suffix)
    except:
        None
    script.close()

    # start simulation
    if socket.gethostname() == 'MAPHost':
        sp.Popen(["bash", scriptname],
                 stdout=open(outname,"w"),stderr=open(errname,"w"))
    else:
        sp.Popen(["sbatch", "-J", jobname, scriptname])

    print "job %s started" %jobname
