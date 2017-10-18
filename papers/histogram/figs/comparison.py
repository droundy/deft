from __future__ import division
import numpy, sys, os
import matplotlib.pyplot as plt
import readnew
from glob import glob

if os.path.exists('../data'):
    os.chdir('..')

energy = int(sys.argv[1])
reference = sys.argv[2]
filebase = sys.argv[3]
methods = [ '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3', '-vanilla_wang_landau']

# For WLTMMC compatibility with LVMC
lvextra = glob('data/%s-wltmmc*-movie' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]

for j in range(len(split2)):
    methods.append('-%s' %split2[j])
print methods

ref = "data/" + reference
maxref = readnew.max_entropy_state(ref)
goodenough = 0.1
#I commented out associated colors for each method as I wasn't sure how
#to modify the colors dict to work for lvmc.  I'll work on fixing this
#in the future.
'''
colors = {
    '-tmi3': 'b',
    '-tmi2': 'k',
    '-tmi': 'y',
    '-toe': 'm',
    '-toe3': 'r',
    '-tmmc': 'k',
    '-toe2': 'c',
    '-vanilla_wang_landau': 'm',
    '-wltmmc': 'c'}
'''
try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)

for method in methods:
    try: 
        r = glob('data/%s%s-movie/*lndos.dat' % (filebase,method))
        if len(r)==0:
            continue
        iterations = numpy.zeros(len(r))
        Nrt_at_energy = numpy.zeros(len(r))
        maxentropystate = numpy.zeros(len(r))
        minimportantenergy = numpy.zeros(len(r))
        erroratenergy = numpy.zeros(len(r))
        goodenoughenergy = numpy.zeros(len(r))
        errorinentropy = numpy.zeros(len(r))
        maxerror = numpy.zeros(len(r))

        for i,f in enumerate(sorted(r)):
            e,lndos,Nrt = readnew.e_lndos_ps(f)

            iterations[i] = readnew.iterations(f)
            Nrt_at_energy[i] = Nrt[energy]
            maxentropystate[i] = readnew.max_entropy_state(f)
            minimportantenergy[i] = readnew.min_important_energy(f)
            doserror = lndos - lndos[maxref] - lndosref + lndosref[maxref]
            errorinentropy[i] = numpy.sum(doserror)/len(doserror)
            erroratenergy[i] = doserror[energy]
            maxerror[i] = numpy.amax(doserror)

        i = 1
        while i < len(iterations) and iterations[i] > iterations[i-1]:
            num_frames_to_count = i+1
            i+=1
        iterations = iterations[:num_frames_to_count]
        minimportantenergy = minimportantenergy[:num_frames_to_count]
        maxentropystate = maxentropystate[:num_frames_to_count]
        Nrt_at_energy = Nrt_at_energy[:num_frames_to_count]
        erroratenergy = erroratenergy[:num_frames_to_count]

        plt.figure('error-at-energy-iterations')
        plt.plot(iterations, erroratenergy, label = method[1:])
        plt.title('Error at energy %g %s' % (energy,filebase))
        plt.xlabel('# iterations')
        plt.ylabel('error')
        plt.legend(loc = 'best')

        plt.figure('round-trips-at-energy' )
        plt.plot(iterations, Nrt_at_energy, label = method[1:])
        plt.title('Roundy Trips at energy %g, %s' % (energy,filebase))
        plt.xlabel('# iterations')
        plt.ylabel('Roundy Trips')
        plt.legend(loc = 'best')
        
        plt.figure('error-at-energy-round-trips')
        plt.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0], label = method[1:])
        plt.title('Error at energy %g %s' % (energy,filebase))
        plt.xlabel('Roundy Trips')
        plt.ylabel('Error')
        plt.legend(loc = 'best')

        plt.figure('maxerror')
        plt.loglog(iterations, maxerror, label = method[1:])
        plt.xlabel('# iterations')
        plt.ylabel('Maximum Entropy Error')
        plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        plt.legend(loc = 'best')

        plt.figure('errorinentropy')
        plt.loglog(iterations, errorinentropy[0:len(iterations)], label = method[1:])
        plt.xlabel('#iterations')
        plt.ylabel('Error in Entropy')
        plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
        plt.legend(loc='best')
    except:
        print 'I had trouble with', method
        raise
plt.show()




















