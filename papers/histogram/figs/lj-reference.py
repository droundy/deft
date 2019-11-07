#!/usr/bin/env python3

from __future__ import division
import sys, os, re
import numpy as np
import readnew
from glob import glob

import yaml

filename = sys.argv[1]

# save the file as
save_name = sys.argv[2]

# Read in energy and lndos
energy = np.loadtxt('%s.energy' % filename)
lndos = np.loadtxt('%s.entropy' % filename)

max_entropy_state = energy[-1]
min_important_energy = energy[0]

print('Smin: ', min_important_energy)
print('Smax: ', max_entropy_state)

dirname = 'data/%s-reference-lndos.dat' % (save_name)
print('saving to', dirname)

np.savetxt(dirname,
          np.c_[energy, lndos[-1]],
          fmt = ('%.16g'),
          delimiter = '\t',
          header = 'comparison reference file\t(generated with python %s \n max_entropy_state: %.16g \n min_important_energy: %.16g \n energy\t lndos\t ' % (' '.join(sys.argv),max_entropy_state,min_important_energy))
