from __future__ import division
import numpy, sys, os, matplotlib

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

methods = [ '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau', '-samc', '-satmmc', '-sad', '-sad3']
# For WLTMMC compatibility with LVMC
lvextra = glob('data/comparison/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]

for j in range(len(split2)):
    methods.append('-%s' %split2[j])
print 'methods are', methods
for method in [mm for m in methods for mm in [m, m+'-tm']]:
    print 'trying method', method
    try:
        dirname = 'data/comparison/%s%s/' % (filebase,method)
        if not os.path.exists(dirname) or os.listdir(dirname) == []:
                continue

        if energy > 0:
                Nrt_at_energy, erroratenergy = numpy.loadtxt(dirname + 'energy-%s.txt' % energy, delimiter = '\t', unpack = True)
        iterations, errorinentropy, maxerror = numpy.loadtxt(dirname + 'errors.txt', delimiter = '\t', unpack = True)

        if not os.path.exists('figs/lv'):
                os.makedirs('figs/lv')
        if energy > 0:
                plt.figure('error-at-energy-iterations')
                colors.plot(iterations, erroratenergy, method=method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('error')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-error-energy-%g.pdf' % (tex_filebase, energy))

                plt.figure('round-trips-at-energy' )
                colors.plot(iterations, Nrt_at_energy, method = method[1:])
                plt.title('Roundy Trips at energy %g, %s' % (energy,filebase))
                plt.xlabel('# iterations')
                plt.ylabel('Roundy Trips')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-round-trips-%g.pdf' % (tex_filebase, energy))

                plt.figure('error-at-energy-round-trips')
                colors.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0],
                         method = method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('Roundy Trips')
                plt.ylabel('Error')
                plt.legend(loc = 'best')
                plt.savefig('figs/%s-error-energy-Nrt-%g.pdf' % (tex_filebase, energy))

        plt.figure('maxerror')
        colors.loglog(iterations, maxerror, method = method[1:])
        plt.xlabel('# iterations')
        plt.ylabel('Maximum Entropy Error')
        plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        plt.legend(loc = 'best')
        plt.savefig('figs/%s-max-entropy-error.pdf' % tex_filebase)

        plt.figure('errorinentropy')
        colors.loglog(iterations, errorinentropy[0:len(iterations)],
                      method = method[1:])
        plt.xlabel('#iterations')
        plt.ylabel('Error in Entropy')
        plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
        plt.legend(loc='best')
        plt.savefig('figs/%s-entropy-error.pdf' % tex_filebase)
    except:
        raise
plt.show()
