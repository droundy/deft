#!/usr/bin/env python2
import matplotlib, sys, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

sys.path.append('../figs')
import styles
import readandcompute

thesis_dir = os.getcwd()
os.chdir('../')

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)


if len(sys.argv) != 6:
    print 'useage: %s ww ff N methods seeds' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = float(sys.argv[3])
#arg N = [25]

methods = eval(sys.argv[4])
#arg methods = [["cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

seeds = eval(sys.argv[5])
#arg seeds = [range(30)]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for seed in seeds for data in ["E","lnw"]]

golden = "tmmc-golden"

max_T = 1.4
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity
S = {} # entropy

# we want to keep our methods distinct from our golden
if golden in methods:
    methods.remove(golden)

fig_size = (6,4.5)
plt.figure('u',figsize=fig_size)
plt.figure('hc',figsize=fig_size)
plt.figure('s',figsize=fig_size)
plt.figure('u_err',figsize=fig_size)
plt.figure('hc_err',figsize=fig_size)
plt.figure('s_err',figsize=fig_size)

array_len = len(readandcompute.u_cv_s(ww, ff, N, methods[0], seeds[0])[0])

for method in [golden]+methods:

    U[method] = numpy.zeros(array_len)
    CV[method] = numpy.zeros(array_len)
    S[method] = numpy.zeros(array_len)
    for seed in seeds:
        u_cv_s = readandcompute.u_cv_s(ww, ff, N, method, seed)
        U[method] += u_cv_s[0] # internal energy
        CV[method] += u_cv_s[1] # heat capacity
        S[method] += u_cv_s[2] # entropy
    U[method] /= len(seeds)
    CV[method] /= len(seeds)
    S[method] /= len(seeds)

    plt.figure('u')
    plt.plot(T_range,U[method]/N, styles.plot(method), label=styles.title(method))

    if method != 'wang_landau':
        plt.figure('hc')
        plt.plot(T_range,CV[method]/N, styles.plot(method), label=styles.title(method))

    plt.figure('s')
    plt.plot(T_range,S[method]/N, styles.plot(method), label=styles.title(method))

methods.remove('wang_landau')

for method in methods:

    plt.figure('u_err')
    plt.plot(T_range,(U[method]-U[golden])/N, styles.plot(method),
             label=styles.title(method))

    plt.figure('hc_err')
    plt.plot(T_range,(CV[method]-CV[golden])/N, styles.plot(method),
             label=styles.title(method))

    plt.figure('s_err')
    plt.plot(T_range,(S[method]-S[golden])/N, styles.plot(method),
             label=styles.title(method))

os.chdir(thesis_dir)

plt.figure('u')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("figs/internal-energy.pdf")

plt.figure('hc')
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("figs/heat-capacity.pdf")

plt.figure('s')
plt.ylim(-7,0)
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textrm{config}}/Nk$')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig("figs/config-entropy.pdf")

plt.figure('u_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='upper right')
plt.ylim(-.01,.02)
plt.tight_layout()
plt.savefig("figs/u-err.pdf")

plt.figure('hc_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta C_V/Nk$')
plt.legend(loc='upper right')
plt.ylim(-0.2,0.2)
plt.tight_layout()
plt.savefig("figs/cv-err.pdf")

plt.figure('s_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel(r'$\Delta S_{\textrm{config}}/Nk$')
plt.legend(loc='upper right')
plt.ylim(-.04,.04)
plt.tight_layout()
plt.savefig("figs/s-err.pdf")
