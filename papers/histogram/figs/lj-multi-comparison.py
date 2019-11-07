#!/usr/bin/env python3

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
#import re

import yaml
import os.path
import time # Need to wait some time if file is being written
np.set_printoptions(threshold=sys.maxsize)

# Example: /home/jordan/sad-monte-carlo/
filename_location = sys.argv[1]

# Example: data/lj-31-reference-lndos.dat
reference = sys.argv[2]

# Used for where we save the data. s000/periodic-ww1.50-ff0.17-N256
filebase = sys.argv[3]

# The number to divide moves by! N is added back in comparison-plot
N = int(sys.argv[4])

seed_avg = int(sys.argv[5])

filename = sys.argv[6:]
print('filenames are ', filename)

# the reference data
ref_data = np.loadtxt(reference,delimiter="\t")
ref_energy = ref_data[:,0]
ref_lndos = ref_data[:,1]

# The Emin and Emax come from the reference energy file rather than
# as command line input.
Emin = ref_energy[-1]
Emax = ref_energy[0]

# We need to deal with float values and indexing integers.
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    # return array[idx] # Just need the index not the value.
    return idx

for f in filename:
    err_in_S = []
    err_max = []
    min_moves = []
    name = '%s.yaml' % (f)
    for n in range(1,seed_avg+1):
        #try:
            name = '%s-s%s' % (f,n)
            print('trying filename ', name)

            # Read the energy and entropy files associated with yaml file.
            if os.path.isfile('%s.yaml' % name):
                energies = np.loadtxt('%s.energy' %  name)
                lndos = np.loadtxt('%s.entropy' %  name)
                moves = np.loadtxt('%s.time' % name)
            else:
                print('There is only a single file with no seeds!\n')
                name = '%s' % (f)
                print('trying new filename ', name)
                if os.path.isfile('%s.yaml' % name):
                    energies = np.loadtxt('%s.energy' %  name)
                    lndos = np.loadtxt('%s.entropy' %  name)
                    moves = np.loadtxt('%s.time' % name)
                    #print(len(energies),len(lndos))
                else:
                    print('unable to read file', name)
                    raise ValueError("\n %s isn't a file!" % ('%s.yaml' % name))

            N_save_times = len(moves)
            # Try to figure out what this was doing and how to reconcile integer and float!
            try:
                maxyaml = Emin
            except:
                my_minE = energies[0]
                num_new_energies = int(my_minE - (-Emin))
                print('num_new_maxyaml', num_new_energies)
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, num_new_energies:] = lndos[:,:]
                lndos = lndos_new
                energies = [0]*num_new_energies + energies
                maxyaml = 0

            try:
                minyaml = Emax
            except:
                my_maxE = energies[-1]
                num_new_energies = -int(my_maxE - (-Emax))
                print('num_new_minyaml', num_new_energies)
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, :lndos.shape[1]] = lndos[:,:]
                lndos = lndos_new
                minyaml = lndos.shape[1]-1

            n_energies = len(ref_energy)
            #print(n_energies)

            errorinentropy = np.zeros(N_save_times)
            maxerror = np.zeros(N_save_times)

            for i in range(0,N_save_times):
                # Compute the norm factor and the error in the DOS.
                yaml_lndos = lndos[i][find_nearest(energies,minyaml+1):find_nearest(energies,maxyaml)]
                yaml_lndos[yaml_lndos == np.inf] = 1e300 # make inf a big number as some early save times have inf!
                yaml_ref_lndos = ref_lndos[find_nearest(energies,minyaml+1):find_nearest(energies,maxyaml)]

                norm_factor = np.mean(yaml_lndos) - np.mean(yaml_ref_lndos)
                doserror = yaml_lndos - yaml_ref_lndos - norm_factor

                errorinentropy[i] = np.sum(abs(doserror))/len(doserror) #- np.mean(doserror)
                maxerror[i] = np.amax(doserror) - np.amin(doserror)


            # remove N from moves in yaml file because N is added back in the
            # comparison-plot script
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

    # for i in range(len(err_in_S)):
    #     err_in_S[i] = err_in_S[i][:len(min_moves)]
    #     err_max[i] = err_max[i][:len(min_moves)]
    # errorinentropy = np.average(err_in_S, axis=0)
    # maxmean = np.amax(err_in_S, axis=0)
    # minmean = np.amin(err_in_S, axis=0)
    # maxerror = np.average(err_max, axis=0)
    # dirname = 'data/comparison/%s-%s' % (filebase, name.replace('-s%s.yaml' %seed_avg,''))
    # print 'saving to', dirname
    # try:
    #     os.mkdir(dirname)
    # except OSError:
    #     pass
    # else:
    #     print ("Successfully created the directory %s " % dirname)
    # np.savetxt('%s/errors.txt' %(dirname),
    #           np.c_[min_moves, errorinentropy, maxerror, minmean, maxmean],
    #           fmt = ('%.4g'),
    #           delimiter = '\t',
    #           header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))
    #
    # # The following is intended for testing whether there is a
    # # systematic error in any of our codes.
    # #np.savetxt('%s/error-vs-energy.txt' %(dirname),
    #             #np.c_[eref, doserror],
    #             #fmt = ('%.4g'),
    #             #delimiter = '\t', header = 'E\t Serror')
