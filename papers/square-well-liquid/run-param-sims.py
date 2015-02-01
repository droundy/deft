#!/usr/bin/python2

import os, sys, socket
import subprocess as sp
import importlib
swmc = importlib.import_module("run-swmc")
import threading

if not len(sys.argv) in [5,6]:
    print 'useage: %s ww ff N method [-to]' % __file__
    exit(1)

method = sys.argv[4]

input_args = sys.argv[1:5]

if '-to' in sys.argv:
    input_args += ["--toggle","['transition_override']"]

threads = []

# build scons before sending off simulations
swdir = os.path.dirname(os.path.realpath(__file__))
projectdir = os.path.realpath(swdir+'/../..')
simname = 'square-well-monte-carlo'

cores = 6 if socket.gethostname() == 'quipu' else 4

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

elif method == 'walker_optimization':
    for init_min_energy_samples in [2,3,5,10]:
        method_args = ["--suffix","params-%i" % init_min_energy_samples,
                       "--values","['init_min_energy_samples','%g']"
                       % init_min_energy_samples]
        threads.append(threading.Thread(target = swmc.run,
                                        args = [input_args + method_args]))
        threads[-1].start()

elif method == 'robustly_optimistic':
    for robust_update_scale in [0.1,0.3,0.5,0.7,0.9]:
        for robust_cutoff in [0.01,0.05,0.1,0.3]:
            method_args = ["--suffix","params-%g-%g" %(robust_update_scale, robust_cutoff),
                           "--values","['robust_update_scale','%g'," % robust_update_scale \
                           + "'robust_cutoff','%g']" % robust_cutoff]
            threads.append(threading.Thread(target = swmc.run,
                                            args = [input_args + method_args]))
            threads[-1].start()


elif method == 'bubble_suppression':
    for bubble_scale in [ 2**n for n in [-1,0,1,2,3] ]:
        for bubble_cutoff in [ 2**n for n in [-3,-2,-1] ]:
            method_args = ["--suffix","params-%g-%g" %(bubble_scale, bubble_cutoff),
                           "--values","['bubble_scale','%g'," % bubble_scale \
                           + "'bubble_cutoff','%g']" % bubble_cutoff]
            threads.append(swmc.run(input_args + method_args))
            threads.append(threading.Thread(target = swmc.run,
                                            args = [input_args + method_args]))
            threads[-1].start()

