#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
#import re

import yaml
import os.path
import time # Need to wait some time if file is being written

# Example: /home/jordan/sad-monte-carlo/
filename_location = sys.argv[1]

# Example: data/samc-1e4-256-cpp-reference-lndos.dat
reference = sys.argv[2]

# Used for where we save the data.: s000/periodic-ww1.50-ff0.17-N256
filebase = sys.argv[3]

# The number to divide moves by! N is added back in comparison-plot
N = int(sys.argv[4])

# Energy range
Emin = int(sys.argv[5])
Emax = int(sys.argv[6])

# Are you comparing to a yaml reference?
yamlRef = bool(sys.argv[7])

seed_avg = int(sys.argv[8])

filename = sys.argv[9:]
print('filenames are ', filename)

for f in filename:
    err_in_S = []
    err_max = []
    min_moves = []
    name = '%s.yaml' % (f)
    for n in range(1,seed_avg+1):
        #try:
            name = '%s-s%s.yaml' % (f,n)
            print('trying filename ', name)

            # Read YAML file
            if os.path.isfile(filename_location + name):
                with open(filename_location + name, 'r') as stream:
                    yaml_data = yaml.load(stream)
            else:
                print('unable to read file', filename_location + name)
                raise ValueError("%s isn't a file!" % (filename_location + name))

            #print(data_loaded)
            data = yaml_data
            data['bins']['histogram'] = np.array(data['bins']['histogram'])
            data['bins']['lnw'] = np.array(data['bins']['lnw'])
            
            data['movies']['entropy'] = np.array(data['movies']['entropy'])
            lndos = data['movies']['entropy']
            N_save_times = len(data['movies']['entropy'])
            try:
                maxyaml = data['movies']['energy'].index(-Emin)
            except:
                my_minE = data['movies']['energy'][0]
                num_new_energies = int(my_minE - (-Emin))
                print('num_new_maxyaml', num_new_energies)
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, num_new_energies:] = lndos[:,:]
                lndos = lndos_new
                maxyaml = 0

            try:
                minyaml = data['movies']['energy'].index(-Emax)
            except:
                my_maxE = data['movies']['energy'][-1]
                num_new_energies = -int(my_maxE - (-Emax))
                print('num_new_minyaml', num_new_energies)
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, :lndos.shape[1]] = lndos[:,:]
                lndos = lndos_new
                minyaml = lndos.shape[1]-1

            #moves = data['moves']
            
            
            ref = reference
            if ref[:len('data/')] != 'data/':
                ref = 'data/' + ref
            maxref = Emax #int(readnew.max_entropy_state(ref))
            minref = Emin # int(readnew.min_important_energy(ref))
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
            moves = data['movies']['time']
            if min_moves == [] or len(min_moves) > len(moves):
                min_moves = np.array(moves)/N
            errorinentropy = errorinentropy[:len(moves)]
            maxerror = maxerror[:len(moves)]
            err_in_S.append(errorinentropy)
            err_max.append(maxerror)

            dirname = 'data/comparison/%s-%s' % (filebase, name.replace('.yaml',''))
            print 'saving to', dirname
            try:  
                os.mkdir(dirname)
            except OSError:  
                pass
            else:  
                print ("Successfully created the directory %s " % dirname)
            np.savetxt('%s/errors.txt' %(dirname),
              np.c_[np.array(moves)/N, errorinentropy, maxerror],
              fmt = ('%.4g'),
              delimiter = '\t',
              header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))
        #except:
        #    print('I ran into some odd trouble.')
        #    pass

    for i in range(len(err_in_S)):
        err_in_S[i] = err_in_S[i][:len(min_moves)]
        err_max[i] = err_max[i][:len(min_moves)]
    errorinentropy = np.average(err_in_S, axis=0)
    maxmean = np.amax(err_in_S, axis=0)
    minmean = np.amin(err_in_S, axis=0)
    maxerror = np.average(err_max, axis=0)
    dirname = 'data/comparison/%s-%s' % (filebase, name.replace('-s%s.yaml' %seed_avg,''))
    print 'saving to', dirname
    try:  
        os.mkdir(dirname)
    except OSError:  
        pass
    else:  
        print ("Successfully created the directory %s " % dirname)
    np.savetxt('%s/errors.txt' %(dirname),
              np.c_[min_moves, errorinentropy, maxerror, minmean, maxmean],
              fmt = ('%.4g'),
              delimiter = '\t',
              header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))
    
    # The following is intended for testing whether there is a
    # systematic error in any of our codes.
    #np.savetxt('%s/error-vs-energy.txt' %(dirname),
                #np.c_[eref, doserror],
                #fmt = ('%.4g'),
                #delimiter = '\t', header = 'E\t Serror')
    
