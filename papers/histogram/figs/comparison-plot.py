from __future__ import division
import sys, os, matplotlib
import numpy as np

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

if os.path.exists('../data'):
    os.chdir('..')

energy = int(sys.argv[1])
filebase = sys.argv[2]
tex_filebase = filebase.replace('.','_') # latex objects to extra "." characters

methods = [ '-sad', '-sad3', '-sad3-s1', '-sad3-s2',
            '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau']
if 'allmethods' not in sys.argv:
    methods = ['-sad3','-tmmc', '-vanilla_wang_landau']
    
# For WLTMMC compatibility with LVMC
lvextra = glob('data/comparison/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]
for meth in split2:
    if meth[-3:] != '-tm':
        methods.append('-%s' % meth)

# For SAMC compatibility with LVMC
lvextra1 = glob('data/comparison/%s-samc*' % filebase)
split3 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra1]
split4 = [i.split('-m', 1)[0] for i in split3]
for meth in split4:
    if meth[-3:] != '-tm':
            methods.append('-%s' % meth)

print 'methods are', methods
for method in methods:
    print 'trying method', method
    try:
        dirname = 'data/comparison/%s%s/' % (filebase,method)
        if not os.path.exists(dirname) or os.listdir(dirname) == []:
                continue

        if energy > 0:
                Nrt_at_energy, erroratenergy = np.loadtxt(dirname + 'energy-%s.txt' % energy, delimiter = '\t', unpack = True)
        iterations, errorinentropy, maxerror = np.loadtxt(dirname + 'errors.txt', delimiter = '\t', unpack = True)

        if not os.path.exists('figs/lv'):
                os.makedirs('figs/lv')
                NxN = filebase.split('-')
                # Formula to calculate N from title i.e. 100x10
                # and use floor to always round up.
                N = np.floor(0.25*0.20*NxN[0]*NxN[-1]*NxN[-1])
                moves = iterations * float(N)

        if not os.path.exists('figs/s000'):
                os.makedirs('figs/s000')
                N = filebase.split('-N')[-1]
                # Get N directly from title.
                moves = iterations * float(N)

        if energy > 0:
                plt.figure('error-at-energy-iterations')
                colors.plot(iterations, erroratenergy, method=method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('error')
                colors.legend()
                plt.savefig('figs/%s-error-energy-%g.pdf' % (tex_filebase, energy))

                plt.figure('round-trips-at-energy' )
                colors.plot(iterations, Nrt_at_energy, method = method[1:])
                plt.title('Round Trips at energy %g, %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('Round Trips')
                colors.legend()
                plt.savefig('figs/%s-round-trips-%g.pdf' % (tex_filebase, energy))

                plt.figure('error-at-energy-round-trips')
                colors.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0],
                         method = method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('Round Trips')
                plt.ylabel('Error')
                colors.legend()
                plt.savefig('figs/%s-error-energy-Nrt-%g.pdf' % (tex_filebase, energy))

        plt.figure('maxerror')
        colors.loglog(iterations, maxerror, method = method[1:])
        plt.xlabel('# iterations')
        plt.ylabel('Maximum Entropy Error')
        plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        colors.legend()
        plt.savefig('figs/%s-max-entropy-error.pdf' % tex_filebase)

        plt.figure('errorinentropy')
        colors.loglog(iterations, errorinentropy[0:len(iterations)],
                      method = method[1:])
        plt.xlabel('#Moves')
        plt.ylabel('Average Entropy Error')
        plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
        colors.legend()
        plt.savefig('figs/%s-entropy-error.pdf' % tex_filebase)

    except:
        raise
plt.show()
