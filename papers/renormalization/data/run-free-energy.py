#!/usr/bin/python
# this should run and plot the free energy simulation eventually

##########
#   Dear next person who has to edit this file,
#   
#   It needs a refactor badly.
#   I'm really sorry.
##########

import numpy as np
import sys, os, math, matplotlib
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino'],'size':'16'})
# rc('text', usetex=True)

# input: "../../free-energy-monte-carlo"
# input: "../../free-energy-monte-carlo-infinite-case"

# cache-suffix: .dat

FILENAME = "periodic-ff%g-ff_small%g-N%d-iterations%d-seed%d-sc_period%d"
data_dir = 'data'

def main(argv=None):
    if argv is None:
        argv = sys.argv

    if 'show' in argv:
        argv.remove('show')

    if len(argv[1:]) == 0:
        run_simulation()
    elif argv[1] == 'triangles':
        check_triangles()
    elif argv[1] == 'dir':
        plot_from_directory(argv[2])
    else:
        print "invalid option"


def run_simulation():
    plot = True
    ffs = []
    success_ratios = []
    all_total_checks = []
    all_valid_checks = []
    steps = 20
    ff = 0
    step_size = 0.05
    N = 10
    sim_iterations = 1000000
    seed = 0
    sc_period = max(10, (1 * N * N)/10)

    for i in xrange(steps):
        # hacky way to reduce step size as we approach a smaller cell
        if i != 0 and i % 4 == 0:
            step_size = step_size * 0.5

        filename = FILENAME %\
            (ff, ff+step_size, N, sim_iterations, seed, sc_period)
        filename_with_extension = filename+"-g.dat"
        filepath = os.path.join(data_dir, filename_with_extension)

        if not os.path.isfile(filepath):
            # then run simulation
            # add args for infinite/regular case, then everything else
            arg_list = []
            if ff == 0:
                arg_list.extend(['../../free-energy-monte-carlo-infinite-case'])
            else:
                arg_list.extend([
                    '../../free-energy-monte-carlo',
                    '--ff', str(ff),
                    '--sc_period', str(sc_period)
                    ])

            arg_list.extend([
                '--iterations', str(sim_iterations),
                '--filename', filename,
                '--data_dir', data_dir,
                #'--seed', str(seed),
                '--ff_small', str(ff+step_size),
                '--N', str(N),
                ])

            print("running with command:")
            print(" ".join(arg_list))
            assert(not subprocess.call(arg_list))

        data = read_data_file_to_dict(filepath)
        next_ff = data['ff_small']
        total_checks = data['total checks of small cell']
        valid_checks = data['total valid small checks']
        success_ratio = (valid_checks * 1.0)/total_checks

        print 'this step had valid_checks', valid_checks
        print 'this step had total_checks', total_checks
        if not ff == 0:
            print 'this step had longest valid run', data['longest valid run']
            print 'this step had longest failed run', data['longest failed run']
            print 'this step had average valid run', data['average valid run']
            print 'this step had average failed run', data['average failed run']
            print 'this step had valid runs', data['valid runs']
            print 'this step had failed runs', data['failed runs']
        print "************ end step ************"

        ffs.append(next_ff)
        success_ratios.append(success_ratio)
        all_total_checks.append(total_checks)
        all_valid_checks.append(valid_checks)

        ff = next_ff

    print ffs
    print success_ratios

    #do plot
    if plot:

        # do error bars
        all_total_checks = np.array(all_total_checks)
        all_valid_checks = np.array(all_valid_checks)

        # ~N =: error in N
        # ~N = sqrt(N)
        #  R = N_a/N_t
        # ~R = ~N_a/N_t
        # Delta F = -kTln(R)
        # ~Delta F = d Delta F/dR ~R
        #          = -kT/R ~N_a/N_t
        #          = -kT/R sqrt(N_a)/N_t = -kt/sqrt(N_a)
        # F = Delta F_1 + Delta F_2 + ... => ~F = sqrt(~F_1^2 + ~F_2^2 + ...)
        error_bars = np.insert(np.sqrt(np.cumsum(1.0/all_valid_checks)), 0, 0)
        print error_bars
        #print map(math.log, error_bars)

        energy = np.insert(-np.cumsum(map(math.log, success_ratios))/N, 0, 0)
        ffs.insert(0, 0)
        ffs = np.array(ffs)
        plt.scatter(ffs, energy, label='Simulation', c='k')
        plt.errorbar(ffs, energy, yerr=error_bars, fmt=None, c='k', ecolor='k')


        plt.plot(ffs, (4*ffs - 3*ffs**2)/(1-ffs)**2, '-', label='Carnahan-Starling')
        plt.xlim(0)
        plt.ylim(0)
        plt.ylabel('$F^C/NkT$')
        plt.xlabel(r'$\eta$')
        plt.legend(loc=2)
        #plt.savefig('rename-me-please.pdf')
        plt.show()


def plot_from_directory(directory = None):
    plot = True
    N_dir = "N10" if directory is None else directory
    N = int(N_dir[1:])
    ffs = []
    success_ratios = []
    all_total_checks = []
    all_valid_checks = []
    P_dir = ""


    for f in filter(lambda f: f.startswith("periodic"), os.listdir(os.path.join(data_dir, N_dir, P_dir))):
        filepath = os.path.join(data_dir, N_dir, P_dir, f)

        data = read_data_file_to_dict(filepath)
        next_ff = data['ff_small']
        total_checks = data['total checks of small cell']
        valid_checks = data['total valid small checks']
        success_ratio = (valid_checks * 1.0)/total_checks

        print 'this step had valid_checks', valid_checks
        print 'this step had total_checks', total_checks
        if ffs:
            print 'this step had longest valid run', data['longest valid run']
            print 'this step had longest failed run', data['longest failed run']
            print 'this step had average valid run', data['average valid run']
            print 'this step had average failed run', data['average failed run']
            print 'this step had valid runs', data['valid runs']
            print 'this step had failed runs', data['failed runs']
        print "************ end step ************"

        ffs.append(next_ff)
        success_ratios.append(success_ratio)
        all_total_checks.append(total_checks)
        all_valid_checks.append(valid_checks)

    print ffs
    print success_ratios

    #do plot
    if plot:

        # do error bars
        all_total_checks = np.array(all_total_checks)
        all_valid_checks = np.array(all_valid_checks)

        # ~N =: error in N
        # ~N = sqrt(N)
        #  R = N_a/N_t
        # ~R = ~N_a/N_t
        # Delta F = -kTln(R)
        # ~Delta F = d Delta F/dR ~R
        #          = -kT/R ~N_a/N_t
        #          = -kT/R sqrt(N_a)/N_t = -kt/sqrt(N_a)
        # F = Delta F_1 + Delta F_2 + ... => ~F = sqrt(~F_1^2 + ~F_2^2 + ...)
        error_bars = np.insert(np.sqrt(np.cumsum(1.0/all_valid_checks)), 0, 0)
        print error_bars
        #print map(math.log, error_bars)

        energy = np.insert(-np.cumsum(map(math.log, success_ratios))/N, 0, 0)
        ffs.insert(0, 0)
        ffs = np.array(ffs)
        plt.scatter(ffs, energy, label='Simulation', c='k')
        plt.errorbar(ffs, energy, yerr=error_bars, fmt=None, c='k', ecolor='k')

        cffs = np.arange(0, ffs[-1], 0.01)
        plt.plot(cffs, (4*cffs - 3*cffs**2)/(1-cffs)**2, '-', label='Carnahan-Starling')
        plt.xlim(0,0.55)
        plt.ylim(0,6)
        plt.ylabel(r'$F^C/NkT$')
        plt.xlabel(r'$\eta$')
        plt.legend(loc=2)
        plt.savefig('cs_%s%s.pdf' % (N_dir, P_dir), bbox_inches='tight')
        #plt.show()


def check_triangles():
    N = 20
    step_sizes = [0.02, 0.01, 0.03]
    ffs = [0.3, 0.32, 0.3]
    sim_iterations = 1000000
    seed = 0
    data_dir = 'data'
    success_ratios = [0] * 3
    sc_period = (1 * N * N)/10


    for i in xrange(3):
        filename = FILENAME %\
            (ffs[i], ffs[i]+step_sizes[i], N, sim_iterations, seed, sc_period)
        filename_with_extension = filename+"-g.dat"
        filepath = os.path.join(data_dir, filename_with_extension)

        if not os.path.isfile(filepath):
            arg_list = [
                '../../free-energy-monte-carlo',
                '--iterations', str(sim_iterations),
                '--filename', filename,
                '--data_dir', data_dir,
                '--seed', str(seed),
                '--ff', str(ffs[i]),
                '--ff_small', str(ffs[i]+step_sizes[i]),
                '--sc_period', str(sc_period)]

            assert(not subprocess.call(arg_list))

        data = read_data_file_to_dict(filepath)
        total_checks = data['total checks of small cell']
        valid_checks = data['total valid small checks']
        success_ratio = (valid_checks * 1.0)/total_checks
        print ffs[i], 'gives valid_checks', valid_checks
        print ffs[i], 'gives total_checks', total_checks
        print ffs[i], 'gives longest valid run', data['longest valid run']
        print ffs[i], 'gives longest failed run', data['longest failed run']
        print ffs[i], 'gives average valid run', data['average valid run']
        print ffs[i], 'gives average failed run', data['average failed run']
        print ffs[i], 'gives valid runs', data['valid runs']
        print ffs[i], 'gives failed runs', data['failed runs']
        success_ratios[i] = success_ratio

    print success_ratios
    print success_ratios[0]*success_ratios[1], 'should be', success_ratios[2]



# I was awake when I wrote this part
def read_data_file_to_dict(filepath):
    print('reading %s\n' % filepath)
    meta_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') and ':' in line:
                sanitized = map(lambda s: s.strip(' #'), line.split(':'))
                meta_data[sanitized[0]] = eval(sanitized[1])
    return meta_data

main()
