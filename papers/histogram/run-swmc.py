#!/usr/bin/python2

import os, sys, socket
import subprocess as sp

def error(args):
    print 'useage: %s ww ff N method --suffix filename_suffix ' % __file__ \
          + '--values value_params --toggle toggle_params'
    print 'args: ', args
    exit(1)

def run(args):

    if len(args) < 4 or len(args) > 10 or len(args)%2 == 1:
        error(args)

    ww = float(args[0])
    ff = float(args[1])
    N = int(args[2])
    method = args[3]

    filename_suffix = ''
    value_params = []
    toggle_params = []

    flag_pos = 4 # position of input flag
    while len(args) > flag_pos:
        if args[flag_pos] == '--suffix':
            filename_suffix = args[flag_pos+1]
        elif args[flag_pos] == '--values':
            value_params = eval(args[flag_pos+1])
        elif args[flag_pos] == '--toggle':
            toggle_params = eval(args[flag_pos+1])
        else:
            error(args)
        flag_pos += 2

    # define some file and directory names
    swdir = os.path.dirname(os.path.realpath(__file__))
    figdir = os.path.realpath(swdir+'/figs')
    projectdir = os.path.realpath(swdir+'/../..')
    jobdir = swdir+'/jobs'
    datadir = swdir+'/data'
    simname = 'square-well-monte-carlo'

    # build SWMC
    os.chdir(swdir)
    exitStatus = sp.call(["fac",simname],
                         stdout = open(os.devnull,"w"),
                         stderr = open(os.devnull,"w"))
    if exitStatus != 0:
        print "fac failed"
        exit(exitStatus)

    memory = 20*N # fixme: better guess
    if 'optimized_ensemble' in toggle_params:
        method += '_oe'
    jobname = 'periodic-ww%04.2f-ff%04.2f-N%i-%s' %(ww, ff, N, method)
    if filename_suffix:
        jobname += '-'+filename_suffix
    basename = "%s/%s" %(jobdir, jobname)
    scriptname = basename + '.sh'
    outname = basename + '.out'
    errname = basename + '.err'

    command = "time nice -19 %s/%s" %(projectdir, simname)

    if not os.path.exists(jobdir):
        os.makedirs(jobdir)
    script = open(scriptname,'w')
    script.write("#!/bin/bash\n")
    script.write("#SBATCH --mem-per-cpu=%i\n" % memory)
    script.write("#SBATCH --output %s\n" % outname)
    script.write("#SBATCH --error %s\n\n" % errname)
    script.write("echo \"Starting job with ID: %s, "
                 "Estimated memory use: %i MB.\"\n\n" %(jobname,memory))
    script.write("cd %s\nnice -19 " %projectdir)
    script.write(command)
    for (arg,val) in [ ("ww",ww), ("ff",ff), ("N",N) ]:
        script.write(" \\\n --%s %s" %(arg,str(val)))
    script.write(" \\\n --%s"%method.replace("kT","kT "))

    if filename_suffix:
        script.write(" \\\n --filename_suffix %s" %filename_suffix)
    for i in range(len(value_params)/2):
        script.write(" \\\n --%s %s" %(value_params[2*i],value_params[2*i+1]))
    for i in range(len(toggle_params)):
        script.write(" \\\n --%s" %toggle_params[i])
    script.write("\n")
    script.close()

    # start simulation
    sp.Popen(["sbatch", "-J", jobname, scriptname])
    print "job %s started" %jobname


if __name__ == "__main__":
    run(sys.argv[1:])
