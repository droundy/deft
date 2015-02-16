#!/usr/bin/python2

import os, sys, socket
import subprocess as sp
import importlib
swmc = importlib.import_module("run-swmc")
import threading

Ns = range(5,31)

if not len(sys.argv) in [5,6]:
    print 'useage: %s ww ff N method' % __file__
    exit(1)

method = sys.argv[4]

input_args = sys.argv[1:]

threads = []

# build scons before sending off simulations
swdir = os.path.dirname(os.path.realpath(__file__))
projectdir = os.path.realpath(swdir+'/../..')
simname = 'square-well-monte-carlo'

cores = 8 if socket.gethostname() == 'quipu' else 4

exitStatus = sp.call(["scons","-j%i"%cores,"-C",projectdir,simname],
                     stdout = open(os.devnull,"w"),
                     stderr = open(os.devnull,"w"))
if exitStatus != 0:
    print "scons failed"
    exit(exitStatus)

if method == 'wang_landau':
    for wl_factor in [ 2**n for n in [-2,-1,0,1] ]:
        for wl_fmod in [1.5, 2, 3]:
            for wl_threshold in [ 2**n for n in [-2,-1,0,1,2] ]:
                for wl_cutoff in [1e-8]:
                    method_args = ["--suffix","params-%g-%g-%g-%g"
                                  %(wl_factor, wl_fmod, wl_threshold, wl_cutoff),
                                  "--values","['wl_factor','%g'," % wl_factor \
                                  + "'wl_fmod','%g'," % wl_fmod \
                                  + "'wl_threshold','%g'," % wl_threshold \
                                  + "'wl_cutoff','%g']" % wl_cutoff]
                    threads.append(threading.Thread(target = swmc.run,
                                                    args = [input_args + method_args]))
                    threads[-1].start()

elif method in ['optimized_ensemble','robustly_optimistic']:
    for init_samples in [5,10,15,20]:
        method_args = ["--suffix","params-%i" % init_samples,
                       "--values","['init_samples','%g']"
                       % init_min_energy_samples]
        threads.append(threading.Thread(target = swmc.run,
                                        args = [input_args + method_args]))
        threads[-1].start()
