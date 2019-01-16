#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
#import re

import yaml

# Example: /home/jordan/sad-monte-carlo/
filename_location = sys.argv[1]

# Example: data/samc-1e4-256-cpp-reference-lndos.dat
reference = sys.argv[2]

# Used for where we save the data.: s000/periodic-ww1.50-ff0.17-N256
filebase = sys.argv[3]

# The number to divide moves by! N is added back in comparison-plot
N = int(sys.argv[4])

# Energy range
Smin = int(sys.argv[5])
Smax = int(sys.argv[6])

# Are you comparing to a yaml reference?
yamlRef = bool(sys.argv[7])

filename = sys.argv[8:]
print('filenames are ', filename)

for f in filename:
    name = '%s.yaml' % (f)
    print('trying filename ', name)
    # Read YAML file
    with open(filename_location + name, 'r') as stream:
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
    #print maxref, minref
    try:
        eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
    except:
        eref, lndosref = readnew.e_lndos(ref)
    
    errorinentropy = np.zeros(N_save_times)
    maxerror = np.zeros(N_save_times)
    
    for i in range(0,N_save_times):
        # below just set average S equal between lndos and lndosref
        if yamlRef:
            # if using yaml as a reference the range is from 0 to len while for C++ the range is 
            # from maxref to minref + 1
            norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[0:(minyaml+1-maxyaml)])
            doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[0:(minyaml+1-maxyaml)] - norm_factor
        else:
            norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[maxref:minref+1])
            doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[maxref:minref+1] - norm_factor
        errorinentropy[i] = np.sum(abs(doserror))/len(doserror)
        maxerror[i] = np.amax(doserror) - np.amin(doserror)
    
    # remove N from moves in yaml file because N is added back in the
    # comparison-plot script
    moves = map(int,data['movies']['time'])
    moves = [x / N for x in moves]
    errorinentropy = errorinentropy[:len(moves)]
    maxerror = maxerror[:len(moves)]
    
    dirname = 'data/comparison/%s-%s' % (filebase, name.replace('.yaml',''))
    print 'saving to', dirname
    try:  
        os.mkdir(dirname)
    except OSError:  
        pass
    else:  
        print ("Successfully created the directory %s " % dirname)
    np.savetxt('%s/errors.txt' %(dirname),
              np.c_[moves, errorinentropy, maxerror],
              fmt = ('%.4g'),
              delimiter = '\t',
              header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))
    
    # The following is intended for testing whether there is a
    # systematic error in any of our codes.
    #np.savetxt('%s/error-vs-energy.txt' %(dirname),
                #np.c_[eref, doserror],
                #fmt = ('%.4g'),
                #delimiter = '\t', header = 'E\t Serror')
    
