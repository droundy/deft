#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, argparse

energy_names = {
    2: 'Phi1',
    3: 'Phi2',
    4: 'Phi3',
    5: 'Ftot',
}
energy_latex = {
    2: r'$\Phi_1$',
    3: r'$\Phi_2$',
    4: r'$\Phi_3$',
    5: r'$F_{tot}$',
}

parser = argparse.ArgumentParser(description='Creates plots.')
parser.add_argument('-w', type=float, action='append')
parser.add_argument('-n', type=float, action='append')
parser.add_argument('-p', action='append')

args = parser.parse_args()
if args.w is None:
    args.w = [0.1]
if args.n is None:
    args.n = [1.0]
if args.p is None:
    args.p = list(energy_names.values())

print(args)

os.system('mkdir -p DATA_LOGS/PLOTS/Plots_n1.0')

def create_fig(dataname, energy_column):
    plt.figure()
    energy_name = energy_names[energy_column]
    print dataname, energy_name
    data = np.loadtxt('DATA_LOGS/%s.dat' % dataname)
    plot_name="DATA_LOGS/PLOTS/%s_vs_%s.png" % (dataname, energy_name)

    dx = data[:,0]
    mcerror = data[:,1]
    energy = data[:,energy_column]

    mce_values = set()
    for e in mcerror:
        mce_values.add(e)
    #print(mce_values)

    for e in sorted(mce_values):
        plt.plot(dx[mcerror==e], energy[mcerror==e], 'x', label='mcerror %s' % e)

    plot_title="%s vs dx for %s" % (energy_name, dataname.replace('_', ' '))
    plt.title(plot_title)
    plt.xlabel(r'$\Delta x$')
    plt.ylabel('%s a.k.a. %s' % (energy_latex[energy_column], energy_name))
    plt.legend()
    #plt.legend(loc='best')

    plt.savefig(plot_name)

def create_plots(dataname):
    for column in [2,3,4,5]:
        if energy_names[column] in args.p:
            create_fig(dataname, column)

for n in args.n:
    for w in args.w:
        create_plots('n%.1f_T2_gw_%.2f_fv_0' % (n, w))

plt.show()
