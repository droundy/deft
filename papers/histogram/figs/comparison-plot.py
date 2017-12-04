from __future__ import division
import numpy, sys, os, matplotlib

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob


if os.path.exists('../data'):
    os.chdir('..')
    
energy = int(sys.argv[1])
filebase = sys.argv[2]

methods = [ '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau', '-samc', '-satmmc', '-sad']
# For WLTMMC compatibility with LVMC
lvextra = glob('data/%s-wltmmc*-movie' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]

for j in range(len(split2)):
    methods.append('-%s' %split2[j])
for method in methods:
    try:
        dirname = 'data/comparison/%s%s/' % (filebase,method)
        if os.listdir(dirname) == []:
            continue
        
        if energy > 0:
                Nrt_at_energy, erroratenergy = numpy.loadtxt(dirname + 'energy-%s.txt' % energy, delimiter = '\t', unpack = True)
        iterations, errorinentropy, maxerror = numpy.loadtxt(dirname + 'errors.txt', delimiter = '\t', unpack = True)
    
        if not os.path.exists('figs/lv'):
                os.makedirs('figs/lv')
        if energy > 0: 
                plt.figure('error-at-energy-iterations')
                plt.plot(iterations, erroratenergy, label = method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('error')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-error-energy-%g.pdf' % (filebase, energy))
    
                plt.figure('round-trips-at-energy' )
                plt.plot(iterations, Nrt_at_energy, label = method[1:])
                plt.title('Roundy Trips at energy %g, %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('Roundy Trips')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-round-trips-%g.pdf' % (filebase, energy))
        
                plt.figure('error-at-energy-round-trips')
                plt.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0], label = method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('Roundy Trips')
                plt.ylabel('Error')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-error-energy-Nrt-%g.pdf' % (filebase, energy))
    
        plt.figure('maxerror')
        plt.loglog(iterations, maxerror, label = method[1:])
        plt.xlabel('# iterations')
        plt.ylabel('Maximum Entropy Error')
        plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        plt.legend(loc = 'best')
        plt.savefig(('figs/%s-max-entropy-error' % (filebase)).replace('.','_')+'.pdf')
        
        plt.figure('errorinentropy')
        plt.loglog(iterations, errorinentropy[0:len(iterations)], label = method[1:])
        plt.xlabel('#iterations')
        plt.ylabel('Error in Entropy')
        plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
        plt.legend(loc='best')
        plt.savefig(('figs/%s-entropy-error' % (filebase)).replace('.','_')+'.pdf')
    except:
        raise
plt.show()
