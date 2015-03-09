#!/usr/bin/python
# this should run and plot the free energy simulation eventually

import numpy as np
import sys, os, math, matplotlib
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from subprocess import call

# input: "../../free-energy-monte-carlo"

def main():
    # run iteration
    # parse meta data
    # start iteration at "shrunk" filling fraction
    # goto start

    # here until line 60 I wrote while half asleep, don't judge me.
    ffs = []
    success_ratios = []

    sf = 0.99
    ff = 0.30

    ffs.append(ff)
    success_ratios.append(1) # for now. should be the absolute ratio

    for i in xrange(20):
        filename = "periodic-ww1.00-ff{0}-N10-sf{1:f}".format(ff, sf)

        os.system('../../free-energy-monte-carlo --data_dir data --sf %g --ff %g --filename %s'
                % (sf, ff, filename))
        # call(['../free-energy-monte-carlo',
        #   '--sf', str(sf),
        #   '--ff', str(ff),
        #   '--filename', filename,
        #   ])


        # at some point I need to chdir so that the files aren't down a directory on accident.
        data = read_data_file_to_dict(os.path.join('data', filename+"-g.dat"))
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
    energy = -np.cumsum(map(math.log, success_ratios))
    print energy
    plt.scatter(ffs, energy)
    ffs = np.array(ffs)
    # FIXME units below?
    plt.plot(ffs, (4*ffs - 3*ffs**2)/(1-ffs)**2, '-')
    plt.ylabel('F/kT')
    plt.xlabel(r'$\eta$')
    plt.savefig('rename-me-please.pdf')
    plt.show()



# I was awake when I wrote this part
def read_data_file_to_dict(filepath):
    meta_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') and ':' in line:
                sanitized = map(lambda s: s.strip(' #'), line.split(':'))
                meta_data[sanitized[0]] = eval(sanitized[1])
    return meta_data

main()
