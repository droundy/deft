#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
#import re

import yaml


filename = sys.argv[1]
reference = sys.argv[2]

#used for where we save the data.
filebase = sys.argv[3]
methodname = sys.argv[4]
N = int(sys.argv[5]) # the number to divide moves by!

#energy range
Smin = int(sys.argv[6])
Smax = int(sys.argv[7])

f = '%s.yaml' % (filename)

# Read YAML file
with open(f, 'r') as stream:
    yaml_data = yaml.load(stream)

#print(data_loaded)
data = yaml_data
data['bins']['histogram'] = np.array(data['bins']['histogram'])
data['bins']['lnw'] = np.array(data['bins']['lnw'])
data['movies']['energy']
minyaml = data['movies']['energy'].index(-Smax)
maxyaml = data['movies']['energy'].index(-Smin)
#print(data['bins']['lnw'])
moves = data['moves']

data['movies']['entropy'] = np.array(data['movies']['entropy'])
lndos = data['movies']['entropy']
N_save_times = len(data['movies']['entropy'])

ref = reference
if ref[:len('data/')] != 'data/':
    ref = 'data/' + ref
maxref = Smax #int(readnew.max_entropy_state(ref))
minref = Smin # int(readnew.min_important_energy(ref))
n_energies = int(minref - maxref+1)
print maxref, minref
try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)

errorinentropy = np.zeros(N_save_times)
maxerror = np.zeros(N_save_times)

for i in range(0,N_save_times):
    # below just set average S equal between lndos and lndosref
    norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[maxref:minref+1])
    doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[maxref:minref+1] - norm_factor
    errorinentropy[i] = np.sum(abs(doserror))/len(doserror)
    maxerror[i] = np.amax(doserror) - np.amin(doserror)

#print(doserror)

moves = map(int,data['movies']['time'])
moves = [x / N for x in moves]
errorinentropy = errorinentropy[:len(moves)]
maxerror = maxerror[:len(moves)]

dirname = 'data/comparison/%s-%s' % (filebase,methodname)
print 'saving to', dirname
try:  
    os.mkdir(dirname)
except OSError:  
    print ("Creation of the directory %s failed" % dirname)
else:  
    print ("Successfully created the directory %s " % dirname)
np.savetxt('%s/errors.txt' %(dirname),
          np.c_[moves, errorinentropy, maxerror],
          fmt = ('%.4g'),
          delimiter = '\t',
          header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))

