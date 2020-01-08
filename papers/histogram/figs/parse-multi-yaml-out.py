#!/usr/bin/env python

from __future__ import division
import sys, os, re
import numpy as np
import readnew
from glob import glob

import yaml
import os.path

# Example: /home/jordan/sad-monte-carlo/
filename_location = sys.argv[1]

N = int(sys.argv[2]) # the number of atoms.

filename = sys.argv[3:]

for f in filename:
    avg_gamma = []
    min_moves = []
    
    name = '%s.yaml' % (f)
    num_seed_files = len(glob('%s/%s-s*.yaml' % (filename_location, f)))
    print((glob('%s/%s-s*.yaml' % (filename_location, f))))
    print(('%s' % num_seed_files))
    for n in range(1, num_seed_files+1):
        try:
            name = '%s-s%s.yaml' % (f, n)
            print(('trying filename ', name))
            # Read YAML file
            if os.path.isfile(filename_location + name):
                with open(filename_location + name, 'r') as stream:
                    data = yaml.load(stream)
            else:
                raise ValueError("%s isn't a file or is being written to!" % (filename_location + name))

            # Parse the method name used.
            lst = data['method']
            method = str(lst).replace("{", "").replace("}", "")
            method = method.replace(':', ',').split(',')[0]
            method = method.replace("'", '')
            print(('the Monte-Carlo method is %s' % method))
            
            if method == 'Sad':
              dirname = 'data/gamma/n%s/%s.dat' % (N, name.replace('.yaml', ''))
              print('saving to', dirname)
              moves = data['movies']['gamma_time']
              gamma = data['movies']['gamma'][:len(moves)]
              avg_gamma.append(gamma)
              np.savetxt(dirname,
                        np.c_[moves, gamma],
                        fmt = ('%.16g'),
                        delimiter = '\t',
                        header = 'comparison reference file\t(generated with python %s \n moves\t gamma\t' % (' '.join(sys.argv)))
            else:
              dirname = 'data/gamma/n%s/%s.txt' % (N, filename)
              print('saving to', dirname)
              moves = np.asarray(data['movies']['gamma_time']).astype(np.float)
              gamma = np.asarray(data['movies']['gamma'][:len(moves)]).astype(np.float)
              avg_gamma.append(gamma)
              np.savetxt(dirname,
                        np.c_[moves, gamma],
                        fmt = ('%.16g'),
                        delimiter = '\t',
                        header = 'comparison reference file\t(generated with python %s \n moves\t gamma\t' % (' '.join(sys.argv)))

            if min_moves == [] or len(min_moves) > len(moves):
                min_moves = np.array(moves)/N
        except:
            pass
    for i in range(len(avg_gamma)):
        avg_gamma[i] = avg_gamma[i][:len(min_moves)]
    gamma = np.average(avg_gamma, axis=0)
    maxmean = np.amax(avg_gamma, axis=0)
    minmean = np.amin(avg_gamma, axis=0)

    dirname = 'data/gamma/n%s/%s.dat' % (N, name.replace('-s%s.yaml' % (num_seed_files), ''))
    print('saving to', dirname)
    np.savetxt(dirname,
               np.c_[min_moves, gamma, maxmean, minmean],
               fmt = ('%.16g'),
               delimiter = '\t',
               header = 'comparison reference file\t(generated with python %s \n min_moves\t gamma\t max_mean_gamma\t min_mean_gamma\t' % (' '.join(sys.argv)))
