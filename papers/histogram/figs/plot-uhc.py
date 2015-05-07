#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles
import readandcompute

if len(sys.argv) != 7:
    print 'useage: %s ww ff N methods golden seed' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = float(sys.argv[3])
#arg N = range(5,31)

methods = eval(sys.argv[4])
#arg methods = [["wang_landau","simple_flat","tmmc","oetmmc"]]

golden = sys.argv[5]
#arg golden = ["tmmc-golden"]

seed = int(sys.argv[6])
#arg seed = [0]

# input: ["data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-tmmc-golden-%s.dat" % (ww, ff, N, data) for data in ["E","lnw"]]

max_T = 1.4
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)
min_T = 0 # we will adjust this

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity
S = {} # entropy

for method in set(methods+[golden]):

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
    plt.plot(T_range,(U[method]-U[golden])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('hc_err')
    plt.plot(T_range,(CV[method]-CV[golden])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('s_err')
    plt.plot(T_range,(S[method]-S[golden])/N,
             styles.plot(method),label=styles.title(method))

plt.figure('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))

plt.figure('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))

plt.figure('s')
plt.title('Configurational entropy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-S.pdf" % (ww*100, ff*100, N))

plt.figure('u_err')
plt.title('Error in specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best')
plt.ylim(-.04,.04)
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u_err.pdf" % (ww*100, ff*100, N))

plt.figure('hc_err')
plt.title('Error in specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta C_V/Nk$')
plt.legend(loc='best')
plt.ylim(-1.2,1.2)
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc_err.pdf" % (ww*100, ff*100, N))

plt.figure('s_err')
plt.title('Error in configurational entropy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.ylim(-.1,.1) # zoom in!
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-s_err.pdf" % (ww*100, ff*100, N))
