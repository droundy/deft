from __future__ import division
import sys, os, matplotlib
import numpy as np

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

if os.path.exists('../data'):
    os.chdir('..')

energy = int(sys.argv[1])
filebase = sys.argv[2]
transcale = sys.argv[3] 

tex_filebase = filebase.replace('.','_') # latex objects to extra "." characters

methods = ['-sad3', '-sad3-s1', '-sad3-s2',
            '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau']
if 'allmethods' not in sys.argv:
    methods = ['-sad3','-tmmc', '-vanilla_wang_landau','-vanilla_wang_landau-minE', '-sad3-test','-sad3-T13','-one_over_t_wang_landau-T13-t',
              '-sad-t2-T13','-sad-t-s1-T13','-sad-t2-s3-T13']
    if transcale == 'slow':
        methods = ['-sad3-slow','-tmmc-slow', '-vanilla_wang_landau-slow','-sad3-T13-slow']
    if transcale == 'fast':
        methods = ['-sad3-fast','-tmmc-fast', '-vanilla_wang_landau-fast']

# For WLTMMC compatibility with LVMC
lvextra = glob('data/comparison/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]

for meth in split2:
    if "default" in transcale:
        if meth[-3:] != '-tm' and "slow" not in meth and "fast" not in meth:
            methods.append('-%s' % meth)
    if "slow" in transcale: 
        if meth[-3:] != '-tm' and "slow" in meth:
            methods.append('-%s' % meth)
    if "fast" in transcale:
        if meth[-3:] != '-tm' and "fast" in meth:
            methods.append('-%s' % meth)

# For SAMC compatibility with LVMC
lvextra1 = glob('data/comparison/%s-samc*' % filebase)
split3 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra1]
split4 = [i.split('-m', 1)[0] for i in split3]

for meth in split4:
    if "default" in transcale:
        if meth[-3:] != '-tm' and "slow" not in meth and "fast" not in meth:
            methods.append('-%s' % meth)
    if "slow" in transcale: 
        if meth[-3:] != '-tm' and "slow" in meth:
            methods.append('-%s' % meth)
    if "fast" in transcale:
        if meth[-3:] != '-tm' and "fast" in meth:
            methods.append('-%s' % meth)


best_ever_max = 1e100
max_time = 0
min_error = 1e200
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
        best_ever_max = min(best_ever_max, maxerror.min())

        if not os.path.exists('figs/lv'):
                os.makedirs('figs/lv')

        if not os.path.exists('figs/s000'):
                os.makedirs('figs/s000')

        if filebase.startswith('lv'):
                NxNsplit = filebase.split('-')
                NxN = NxNsplit[-1].split('x')
                # Formula to calculate N from title i.e. 100x10
                # and use floor to always round up.
                N = np.floor(0.25*0.20*float(NxN[0])*float(NxN[-1])*float(NxN[-1]))
                ff = 1.0 # FIXME
                moves = iterations * N
        if filebase.startswith('s000'):
                N = filebase.split('-N')[-1]
                ff = filebase.split('-ff')[-1].split('-N')[0]
                ff = float(ff)
                # Get N directly from title.
                moves = iterations * float(N)
        max_time = max(max_time, moves.max())

        if energy > 0:
                plt.figure('error-at-energy-iterations')
                colors.plot(moves, erroratenergy, method=method[1:])
                plt.title('Error at energy %g %s' % (energy,filebase))
                plt.xlabel('# moves')
                plt.ylabel('error')
                colors.legend()
                plt.savefig('figs/%s-error-energy-%g.pdf' % (tex_filebase, energy))

                plt.figure('round-trips-at-energy' )
                colors.plot(moves, Nrt_at_energy, method = method[1:])
                plt.title('Round Trips at energy %g, %s' % (energy,filebase))
                plt.xlabel('# moves')
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
        colors.loglog(moves, maxerror, method = method[1:])
        plt.xlabel(r'Moves')
        plt.ylabel(r'Maximum Entropy Error')
        #plt.title('Maximum Entropy Error vs Iterations, %s' %filebase)
        colors.legend()

        if type(iterations) is not np.float64:
                plt.figure('errorinentropy')
                my_S_error = errorinentropy[0:len(iterations)]
                min_error = min(min_error, my_S_error[my_S_error > 0].min())
                colors.loglog(moves, my_S_error, method = method[1:])
                plt.xlabel(r'$\textrm{Moves}$')
                plt.ylabel(r'$\textrm{Average Entropy Error}$')
                #plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
                if "default" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.05$' % (int(N), ff))
                elif "slow" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.005$' % (int(N), ff))
                elif "fast" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.5$' % (int(N), ff))
                colors.legend()

    except:
        raise
plt.figure('maxerror')
colors.loglog([1e6,max_time], [best_ever_max*np.sqrt(max_time/1e6), best_ever_max], method = r'$\frac{1}{\sqrt{t}}$')
colors.legend()
plt.savefig('figs/%s-max-entropy-error-%s.pdf' % (tex_filebase,transcale))
 
plt.figure('errorinentropy')
moves = np.array([1e5, 1e13])
#colors.loglog(moves, min_error*np.sqrt(moves.max())/np.sqrt(moves), method = r'1/sqrt(t)')
for i in np.arange(-8, 9, 1.0):
    colors.loglog(moves, 10**i*np.sqrt(moves.max())/np.sqrt(moves), method = r'1/sqrt(t)')
plt.xlim(moves[0], moves[1])
if filebase == 's000/periodic-ww1.30-ff0.30-N50':
    plt.ylim(1e-3, 1e5)
elif filebase == 's000/periodic-ww1.30-ff0.30-N500':
    plt.ylim(1e-2, 1e4)
colors.legend()
plt.tight_layout()
plt.savefig('figs/%s-entropy-error-%s.pdf' % (tex_filebase,transcale))


plt.show()
