#!/usr/bin/env python

from __future__ import division
import sys, os, re
import numpy as np
import readnew
from glob import glob

import yaml

filename = sys.argv[1]

#used for where we save the data.
filebase = sys.argv[2]

#energy range
Smin = int(sys.argv[3])
Smax = int(sys.argv[4])

f = '%s.yaml' % (filename)
filename = filename.split('/')[-1]

# Read YAML file.
with open(f, 'r') as stream:
    yaml_data = yaml.load(stream)

#print(data_loaded)
data = yaml_data

#parse the method name used.
lst = data['method']
method = str(lst).replace("{","").replace("}", "")
method = method.replace(':', ',').split(',')[0]
method = method.replace("'", '')
print 'the Monte-Carlo method is ', method

data['movies']['energy']
minE = data['movies']['energy'].index(-Smax)
maxE = data['movies']['energy'].index(-Smin)
#print minE, maxE
moves = data['movies']['time']

dirname = 'data/%s-reference-lndos.dat' % (filename)
print 'saving to', dirname

well_width = data['system']['cell']['well_width']
translation_scale = ['translation_scale']
if method == 'Sad':
    min_T = data['method'][method]['min_T']
    too_lo = data['method'][method]['too_lo']
    too_hi = data['method'][method]['too_hi']
moves = data['moves']
x = data['system']['cell']['box_diagonal']['x']
y = data['system']['cell']['box_diagonal']['y']
z = data['system']['cell']['box_diagonal']['z']
gamma = data['movies']['gamma'][-1]

# shouldn't need to restric range [maxE:minE+1][::-1]
energy = data['movies']['energy'][maxE:minE+1][::-1]
energy[:] = [x*(-1) for x in energy]
lndos = data['movies']['entropy'][-1][maxE:minE+1][::-1]
lndos[:] = [x - max(lndos) for x in lndos]
#print data['movies']['entropy']
ps = data['bins']['round_trips'][maxE:minE+1][::-1]
#print energy
max_entropy_state = energy[0]
min_important_energy = energy[-1]
np.savetxt(dirname,
          np.c_[energy, lndos, ps],
          fmt = ('%.16g'),
          delimiter = '\t',
          #newline = '# seed:',
          #newline = '# well_width: %g' % (well_width),
          #newline = '# ff: ',
          #newline = '# N: ',
          #newline = '# walls: ',
          #newline = '# cell dimensions: (%g, %g, %g)' % (x,y,z),
          #newline = '# translation_scale: %g' % (translation_scale),
          #newline = '# energy_levels: ',
          #newline = '# min_T: %g' % (min_T),
          #newline = '# max_entropy_state: ',
          #newline = '# min_important_energy: ',
          #newline = '# too_high_energy: %i' % (too_hi),
          #newline = '# too_low_energy: %i' % (too_lo),
          #newline = '',
          ##newline = '# WL Factor: %g' % (gamma),
          ##newline = '# iterations: %i' % (),
          ##newline = '# working moves: %i' % (),
          ##newline = '# total moves: %i' % (),
          ##newline = '# acceptance rate: %i' % (),
          #newline = '',
          ##newline = '# converged state: %i' % (),
          #newline = '# converged temperature: ',
          #newline = '# energy\t lndos\t ps\t lndos_tm: ',
          #newline = '\n version: created with yaml\n',
          header = 'comparison reference file\t(generated with python %s \n max_entropy_state: %i \n min_important_energy: %i \n energy\t lndos\t\t ps\t ' % (' '.join(sys.argv),max_entropy_state,min_important_energy))
