#!/usr/bin/python
# this should run and plot the free energy simulation eventually

import numpy as np
import sys, os, math, matplotlib
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess

# input: "../../free-energy-monte-carlo"
# input: "../../free-energy-monte-carlo-infinite-case"

# cache-suffix: .dat

FILENAME = "periodic-ff%g-ff_small%g-N%d-iterations%d-seed%d"

def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv[1:]) == 0:
        run_simulation()
    elif argv[1] == 'triangles':
        check_triangles()
    else:
        print "invalid option"


def run_simulation():
    ffs = []
    success_ratios = []
    steps = 20
    ff = 0
    step_size = 0.05
    N = 10
    sim_iterations = 1000000
    seed = 0
    data_dir = 'data'

    for i in xrange(steps):
        # hacky way to reduce step size as we approach a smaller cell
        if i != 0 and i % 4 == 0:
            step_size = step_size * 0.5

        filename = FILENAME %\
            (ff, ff+step_size, N, sim_iterations, seed)
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
                    '--ff', str(ff)
                    ])

            arg_list.extend([
                '--iterations', str(sim_iterations),
                '--filename', filename,
                '--data_dir', data_dir,
                #'--seed', str(seed),
                '--ff_small', str(ff+step_size)
                ])

            subprocess.call(arg_list)

        data = read_data_file_to_dict(filepath)
        next_ff = data['ff_small']
        total_checks = data['total checks of small cell']
        valid_checks = data['total valid small checks']
        success_ratio = (valid_checks * 1.0)/total_checks

        ffs.append(next_ff)
        success_ratios.append(success_ratio)

        ff = next_ff

    print ffs
    print success_ratios

    #do plot
    energy = -np.cumsum(map(math.log, success_ratios))/N
    print energy
    plt.scatter(ffs, energy)
    ffs = np.array(ffs)
    # FIXME units below?
    plt.plot(ffs, (4*ffs - 3*ffs**2)/(1-ffs)**2, '-')
    plt.ylabel('F/NkT')
    plt.xlabel(r'$\eta$')
    plt.savefig('rename-me-please.pdf')
    plt.show()


def check_triangles():
    print "checking that two shorter simulations reach the same conclusion as one longer one"

    N = 10
    step_size = 0.025
    ffs = []
    sim_iterations = 1000000
    seed = 0
    data_dir = 'data'

    # filling fractions as vertices
    v1 = .2
    v2 = .25
    v3 = .30

    legs = [0] * 3

    for i, (start, end) in enumerate([(v1,v2), (v2,v3), (v1,v3)]):
        success_ratios = []
        ff = start
        while (ff < end):
            filename = FILENAME %\
                (ff, ff+step_size, N, sim_iterations, seed)
            filename_with_extension = filename+"-g.dat"
            filepath = os.path.join(data_dir, filename_with_extension)

            if not os.path.isfile(filepath):
                arg_list = [
                    '../../free-energy-monte-carlo',
                    '--iterations', str(sim_iterations),
                    '--filename', filename,
                    '--data_dir', data_dir,
                    #'--seed', str(seed),
                    '--ff_small', str(ff+step_size)
                    ]

                subprocess.call(arg_list)

            data = read_data_file_to_dict(filepath)
            next_ff = data['ff_small']
            total_checks = data['total checks of small cell']
            valid_checks = data['total valid small checks']
            success_ratio = (valid_checks * 1.0)/total_checks

            success_ratios.append(success_ratio)
            ff = next_ff
        print success_ratios
        legs[i] =  -np.cumsum(map(math.log, success_ratios))/N

    print legs
    print legs[0][-1] + legs[1][-1]
    print legs[2][-1]


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
