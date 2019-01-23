#!/usr/bin/env python

from __future__ import division
import sys, os, re
import numpy as np
import readnew
from glob import glob

import yaml

filename = sys.argv[1]
N = sys.argv[2] # the number of atoms.

f = '%s.yaml' % (filename)
filename = filename.split('/')[-1]

# Read YAML file.
with open(f, 'r') as stream:
    yaml_data = yaml.load(stream)
data = yaml_data

#parse the method name used.
lst = data['method']
method = str(lst).replace("{","").replace("}", "")
method = method.replace(':', ',').split(',')[0]
method = method.replace("'", '')
print 'the Monte-Carlo method is ', method

if method == 'Sad':
    dirname = 'data/gamma/n%s/%s.dat' % (N,filename)
    print 'saving to', dirname
    #min_T = data['method'][method]['min_T'] #gamma plot requires this as an arg.
    moves = data['movies']['gamma_time']
    #too_hi = data['method'][method]['too_hi']
    #too_lo = data['method'][method]['too_lo']
    gamma = data['movies']['gamma']
    np.savetxt(dirname,
          np.c_[moves,gamma],
          fmt = ('%.16g'),
          delimiter = '\t',
          header = 'comparison reference file\t(generated with python %s \n moves\t gamma\t' % (' '.join(sys.argv)))
else:
    dirname = 'data/gamma/n%s/%s.txt' % (N,filename)
    print 'saving to', dirname
    moves = np.asarray(data['movies']['gamma_time']).astype(np.float)
    gamma = np.asarray(data['movies']['gamma']).astype(np.float)
    np.savetxt(dirname,
          np.c_[moves,gamma],
          fmt = ('%.16g'),
          delimiter = '\t',
          header = 'comparison reference file\t(generated with python %s \n moves\t gamma\t' % (' '.join(sys.argv)))
#print "I don't know what method are you using?"
