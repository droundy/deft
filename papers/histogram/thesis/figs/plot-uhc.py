#!/usr/bin/python2
import matplotlib, sys, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

sys.path.append('../figs')
import styles
import readandcompute

thesis_dir = os.getcwd()
os.chdir('../')

fig_size = numpy.array([5,4])

if len(sys.argv) != 6:
    print 'useage: %s ww ff N methods seed' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = float(sys.argv[3])
#arg N = [25]

methods = eval(sys.argv[4])
#arg methods = [["tmmc-golden","cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

seed = int(sys.argv[5])
#arg seed = [0]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

reference = "tmmc-golden"

max_T = 1.4
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)
min_T = 0 # we will adjust this

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity
S = {} # entropy

# we want to keep our methods distinct from our reference
if reference in methods:
    methods.remove(reference)

plt.figure('u',figsize=fig_size)
plt.figure('hc',figsize=fig_size)
plt.figure('s',figsize=fig_size)
plt.figure('u_err',figsize=fig_size)
plt.figure('hc_err',figsize=fig_size)
plt.figure('s_err',figsize=fig_size)

for method in [reference]+methods:

    u_cv_s = readandcompute.u_cv_s(ww, ff, N, method)
    if u_cv_s != None:
        U[method] = u_cv_s[0] # internal energy
        CV[method] = u_cv_s[1] # heat capacity
        S[method] = u_cv_s[2] # entropy

    plt.figure('u')
    plt.plot(T_range,U[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('hc')
    plt.plot(T_range,CV[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('s')
    plt.plot(T_range,S[method]/N,styles.plot(method),label=styles.title(method))

for method in methods:

    plt.figure('u_err')
    plt.plot(T_range,(U[method]-U[reference])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('hc_err')
    plt.plot(T_range,(CV[method]-CV[reference])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('s_err')
    plt.plot(T_range,(S[method]-S[reference])/N,
             styles.plot(method),label=styles.title(method))

os.chdir(thesis_dir)

plt.figure('u')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/internal-energy.pdf")

plt.figure('hc')
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/heat-capacity.pdf")

plt.figure('s')
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/config-entropy.pdf")

plt.figure('u_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best')
plt.ylim(-.04,.04)
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/internal-energy-err.pdf")

plt.figure('hc_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta C_V/Nk$')
plt.legend(loc='best')
plt.ylim(-1.2,1.2)
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/heat-capacity-err.pdf")

plt.figure('s_err')
plt.xlabel('$kT/\epsilon$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.ylim(-.1,.1) # zoom in!
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/config-entropy-err.pdf")
