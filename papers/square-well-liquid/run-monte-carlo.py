#!/usr/bin/python2

import os, sys, socket
import subprocess as sp

if not len(sys.argv) in [5,6,7]:
    print 'useage: %s ww ff N methods filename_suffix optional_params' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])

ff = float(sys.argv[2])

N = int(sys.argv[3])

methods = eval(sys.argv[4])

if len(sys.argv) >= 6:
    filename_suffix = sys.argv[5]
else:
    filename_suffix = ''

if len(sys.argv) == 7:
    optional_params = eval(sys.argv[6])
else:
    optional_params = []

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

for method in methods:
    memory = 20*N # fixme: better guess
    jobname = 'periodic-ww%04.2f-ff%04.2f-N%i-%s' %(ww, ff, N, method)
    if filename_suffix:
        jobname += '-'+filename_suffix
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
    for (arg,val) in [ ("ww",ww), ("ff",ff), ("N",N) ]:
        script.write(" \\\n --%s %s" %(arg,str(val)))
    script.write(" \\\n --%s"%method.replace("kT","kT "))

    if filename_suffix:
        script.write(" \\\n --filename_suffix %s" %filename_suffix)
    for i in range(len(optional_params)/2):
        script.write(" \\\n --%s %s" %(optional_params[2*i],optional_params[2*i+1]))
    script.close()

    # start simulation
    if socket.gethostname() == 'MAPHost':
        sp.Popen(["bash", scriptname],
                 stdout = open(outname,"w"), stderr = open(errname,"w"))
    else:
        sp.Popen(["sbatch", "-J", jobname, scriptname])

    print "job %s started" %jobname
