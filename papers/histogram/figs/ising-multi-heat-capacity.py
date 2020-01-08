#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import readnew
from glob import glob
#import re
import matplotlib.pyplot as plt
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
print(('filenames are ', filename))

for f in filename:
    err_in_S = []
    err_max = []
    min_moves = []
    name = '%s.yaml' % (f)
    for n in range(1, seed_avg+1):
        #try:
            name = '%s-s%s.yaml' % (f, n)
            print(('trying filename ', name))

            # Read YAML file
            if os.path.isfile(filename_location + name):
                with open(filename_location + name, 'r') as stream:
                    yaml_data = yaml.load(stream)
            else:
                print(('unable to read file', filename_location + name))
                raise ValueError("%s isn't a file!" % (filename_location + name))

            data = yaml_data
            data['bins']['histogram'] = np.array(data['bins']['histogram'])
            data['bins']['lnw'] = np.array(data['bins']['lnw'])

            data['movies']['entropy'] = np.array(data['movies']['entropy'])
            lndos = data['movies']['entropy']
            energies = data['movies']['energy']
            N_save_times = len(data['movies']['entropy'])
            try:
                maxyaml = energies.index(-Emin)
            except:
                my_minE = energies[0]
                num_new_energies = int(my_minE - (-Emin))
                print(('num_new_maxyaml', num_new_energies))
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, num_new_energies:] = lndos[:,:]
                lndos = lndos_new
                energies = [0]*num_new_energies + energies
                maxyaml = 0

            try:
                minyaml = energies.index(-Emax)
            except:
                my_maxE = energies[-1]
                num_new_energies = -int(my_maxE - (-Emax))
                print(('num_new_minyaml', num_new_energies))
                lndos_new = np.zeros((lndos.shape[0], lndos.shape[1]+num_new_energies))
                lndos_new[:, :lndos.shape[1]] = lndos[:,:]
                lndos = lndos_new
                minyaml = lndos.shape[1]-1

            #moves = data['moves']

            # Do the reference calculations:
            ref = reference
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
            # Internal energy and heat capacity error:
            U = np.zeros(N_save_times)
            CV = np.zeros(N_save_times)
            beta = np.linspace(0.2, 1.0, N_save_times)

            for i in range(0, N_save_times):
                # below just set average S equal between lndos and lndosref
                if yamlRef:
                    # if using yaml as a reference the range is from 0 to len while for C++ the range is
                    # from maxref to minref + 1
                    if 'ising' in filebase:

                        ising_norm = lndos[i][maxyaml:minyaml+1] # remove impossible state

                        ising_lndos = lndos[i][maxyaml:minyaml+1][::-1] # remove impossible state

                        # the states are counted backward hence the second to last state would be at index = 1
                        ising_norm = np.delete(ising_norm, [1])
                        ising_lndos = np.delete(ising_lndos, [len(ising_lndos)-2])

                        # Calculate the internal energy and heat capacity:
                        PartitionF = []
                        U_numer = []
                        C_numer = []

                        # Compute the Boltzmann factor and reference
                        Boltz_fac = []
                        Boltz_fac_ref = []

                        ising_E = np.array(energies[maxyaml:minyaml+1])
                        ising_E = np.delete(ising_E, [len(ising_E)-2])
                        for EE in range(len(ising_E)):
                            #print(ising_lndos,ising_lndos.size,ising_E.size,beta[i])
                            Boltz_fac.append(ising_lndos[EE] - beta[i]*ising_E[EE])
                            Boltz_fac_ref.append(lndosref[EE] - beta[i]*ising_E[EE])

                        # Compute the Maximum in the Boltz_fac
                        Max_Boltz_fac = max(Boltz_fac)
                        Max_Boltz_fac_ref = max(Boltz_fac_ref)
                        #print(Max_Boltz_fac, Max_Boltz_fac_ref)

                        for EE in range(len(ising_E)):
                            # print(len(new_ising_lndos))
                            # print(len(ising_E))
                            PartitionF.append(np.array(np.exp(Boltz_fac[EE]-Max_Boltz_fac)))
                        for EE in range(len(ising_E)):
                            U_numer.append(np.array(ising_E[EE]*np.exp(Boltz_fac[EE]-Max_Boltz_fac)))
                        for EE in range(len(ising_E)):
                            C_numer.append(np.array((ising_E[EE]**2)*np.exp(Boltz_fac[EE]-Max_Boltz_fac)))

                        PartitionF = np.asarray(PartitionF)
                        U_numer = np.asarray(U_numer)
                        C_numer = np.asarray(C_numer)
                        #
                        norm_factor = np.mean(ising_norm) - np.mean(lndosref[0:minref-maxref+1])
                        doserror = ising_lndos - lndosref[0:minref-maxref+1] - norm_factor
                    else:
                        norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[0:minref-maxref+1])
                        doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[0:minref-maxref+1] - norm_factor
                else:
                    norm_factor = np.mean(lndos[i][maxyaml:minyaml+1]) - np.mean(lndosref[maxref:minref+1])
                    doserror = lndos[i][maxyaml:minyaml+1][::-1] - lndosref[maxref:minref+1] - norm_factor

                errorinentropy[i] = np.sum(abs(doserror))/len(doserror) #- np.mean(doserror)
                maxerror[i] = np.amax(doserror) - np.amin(doserror)

                U[i] = np.sum(U_numer) / np.sum(PartitionF)
                CV[i] = beta[i]**2 * (np.sum(C_numer) / np.sum(PartitionF) - U[i]**2)

                print(('The error in the internal energy is: ', U[i]))
                print(('The error in the heat capacity is: ', CV[i]))
            plt.figure()
            plt.plot(beta, CV)
            #plt.plot(beta, U)
            plt.show()
            # remove N from moves in yaml file because N is added back in the
            # comparison-plot script
            moves = data['movies']['time']
            if min_moves == [] or len(min_moves) > len(moves):
                min_moves = np.array(moves)/N
            errorinentropy = errorinentropy[:len(moves)]
            maxerror = maxerror[:len(moves)]
            err_in_S.append(errorinentropy)
            err_max.append(maxerror)

            # dirname = 'data/comparison/%s-%s' % (filebase, name.replace('.yaml',''))
            # print 'saving to', dirname
            # try:
            #     os.mkdir(dirname)
            # except OSError:
            #     pass
            # else:
            #     print ("Successfully created the directory %s " % dirname)
            # np.savetxt('%s/errors.txt' %(dirname),
            #   np.c_[np.array(moves)/N, errorinentropy, maxerror],
            #   fmt = ('%.4g'),
            #   delimiter = '\t',
            #   header = 'iterations\t errorinentropy\t maxerror\t(generated with python %s' % ' '.join(sys.argv))
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
    dirname = 'data/comparison/%s-%s' % (filebase, name.replace('-s%s.yaml' %seed_avg, ''))
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
