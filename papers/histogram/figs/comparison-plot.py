from __future__ import division
import sys, os, matplotlib
import numpy as np

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

from matplotlib.colors import LightSource

densitycolormap = plt.cm.jet
densityinterpolation = 'bilinear'
densityshadedflag = True
densitybarflag = True
gridflag = True

if os.path.exists('../data'):
    os.chdir('..')

energy = int(sys.argv[1])
filebase = sys.argv[2]
tex_filebase = filebase.replace('.','_') # latex objects to extra "." characters

methods = [ '-sad', '-sad3', '-sad3-s1', '-sad3-s2',
            '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau']
# For WLTMMC compatibility with LVMC
lvextra = glob('data/comparison/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]
for j in range(len(split2)):
    methods.append('-%s' %split2[j])

# For SAMC compatibility with LVMC
lvextra1 = glob('data/comparison/%s-samc*' % filebase)
split3 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra1]
split4 = [i.split('-m', 1)[0] for i in split3]
for j in range(len(split4)):
    methods.append('-%s' %split4[j])

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

        #if os.path.isfile(dirname + 'wl-factor.txt'):
            #iterations, wl_factor = np.loadtxt(dirname + 'wl-factor.txt', delimiter = '\t', unpack = True)


        if not os.path.exists('figs/lv'):
                os.makedirs('figs/lv')

        if not os.path.exists('figs/s000'):
                os.makedirs('figs/s000')
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

                #------------------------------------------#
                # Perhaps make a subplot for each method? For now I will test on
                # just one Monte-Carlo Method: TMI3.

                #if method == methods[3]:
                        #plt.figure('Maxerror-at-energy-round-trips')
                        #plt.title('Error at energy %g %s' % (energy,filebase))
                        #plt.xlabel('Round Trips')
                        #plt.ylabel('iterations')

                        #X,Y = np.meshgrid(Nrt_at_energy[Nrt_at_energy > 0], iterations[Nrt_at_energy > 0])
                        #Z = np.log(X) + np.log(Y) # Z in no way coresponds to MaxError!
                        #if densityshadedflag:
                                #ls = LightSource(azdeg=120,altdeg=65)
                                #rgb = ls.shade(Z,densitycolormap)
                                #im = plt.imshow(rgb, vmin = Z.min()/2, vmax = 2*Z.max(), cmap=densitycolormap)
                                #cset = plt.contour(Z, np.arange(Z.min(), Z.max(), (Z.max()-Z.min())/6),
                                                   #linewidths=2, cmap=plt.cm.Set2)
                                #plt.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
                        #if densitybarflag:
                                #plt.colorbar(im)
                        #plt.savefig('figs/%s-Maxerror-energy-Nrt-%g.pdf' % (tex_filebase, energy))

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
        plt.xlabel('#iterations')
        plt.ylabel('Error in Entropy')
        plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
        colors.legend()
        plt.savefig('figs/%s-entropy-error.pdf' % tex_filebase)

        #if os.path.isfile(dirname + 'wl-factor.txt'):
            #plt.figure('wl-factor')
            #colors.loglog(iterations, wl_factor,
                      #method = method[1:])
            #plt.xlabel('#iterations')
            #plt.ylabel('WL factor')
            #plt.title('WL Factor at Each Iteration, %s' %filebase)
            #colors.legend()
            #plt.savefig('figs/%s-wl-factor.pdf' % tex_filebase)
    except:
        raise
plt.show()
