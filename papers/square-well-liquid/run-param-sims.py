#!/usr/bin/python2

import os, sys, socket
import subprocess as sp

if not len(sys.argv) in [5,6]:
    print 'useage: %s ww ff N method transition_override' % sys.argv[0]
    exit(1)

method = sys.argv[4]

transition_override = len(sys.argv) == 6

sim_options = ["./run-monte-carlo.py",
               sys.argv[1], sys.argv[2], sys.argv[3], "%s" %method]

exec_list = []

if method == 'wang_landau':
    for wl_factor in [ 2**n for n in [-2,-1,0,1] ]:
        for wl_fmod in [1.5, 2, 3]:
            for wl_threshold in [ 2**n for n in [-2,-1,0,1,2] ]:
                for wl_cutoff in [1e-8]:
                    exec_list.append(
                        sim_options + ["params-%g-%g-%g-%g"
                                       %(wl_factor, wl_fmod, wl_threshold, wl_cutoff),
                                       "['wl_factor','%g'," %wl_factor \
                                       + "'wl_fmod','%g'," %wl_fmod \
                                       + "'wl_threshold','%g',"%wl_threshold \
                                       + "'wl_cutoff','%g']" %wl_cutoff])

elif method == 'walker_optimization':
    for init_min_energy_samples in [5,10,15]:
        exec_list.append(
            sim_options + ["params-%i" %init_min_energy_samples,
                           "['init_min_energy_samples','%g']"
                           %init_min_energy_samples])

elif method == 'robustly_optimistic':
    for robust_update_scale in [0.1,0.3,0.5,0.7,0.9]:
        for robust_cutoff in [0.01,0.05,0.1,0.3]:
            exec_list.append(
                sim_options + ["params-%g-%g" %(robust_update_scale, robust_cutoff),
                               "['robust_update_scale','%g'," %robust_update_scale \
                               + "'robust_cutoff','%g']" %robust_cutoff])

elif method == 'bubble_suppression':
    for bubble_scale in [ 2**n for n in [-1,0,1,2,3] ]:
        for bubble_cutoff in [ 2**n for n in [-3,-2,-1] ]:
            exec_list.append(
                sim_options + ["params-%g-%g" %(bubble_scale, bubble_cutoff),
                               "['bubble_scale','%g'," %bubble_scale \
                               + "'bubble_cutoff','%g']" %bubble_cutoff])

for i in range(len(exec_list)):
    if transition_override:
        exec_list[i] += ['["transition_override"]']
    sp.Popen(exec_list[i])
