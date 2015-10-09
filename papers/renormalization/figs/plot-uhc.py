#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles
import readandcompute

if len(sys.argv) != 5:
    print 'useage: %s i ww L N ' % sys.argv[0]
    exit(1)

i = float(sys.argv[1])

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[3])
#arg L = [5]

N = float(sys.argv[4])
#arg N = range(5,31)

#methods = eval(sys.argv[4])
methods = ["tmmc", "oetmmc"]
#arg methods = [["tmmc"]]

#seed = int(sys.argv[5])
#arg seed = [0]

# input: ["../data/ww%04.2f-L%04.2f-N%i-%s.dat" % (ww, L, N, data) for data in ["E","lnw"]]

for method in set(methods):

    T, U, CV, S, min_T = readandcompute.T_u_cv_s_minT('../data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d' % (ww, L,  i, N))

    plt.figure('u')
    plt.plot(T,U/N,styles.plot(method),label=styles.title(method))

    plt.figure('hc')
    plt.plot(T,CV/N,styles.plot(method),label=styles.title(method))

    plt.figure('s')
    plt.plot(T,S/N,styles.plot(method),label=styles.title(method))

plt.figure('u')
plt.title('Specific internal energy for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("ww%02.0f-L%02.0f-N%i-u.pdf" % (ww*100, L, N))

plt.figure('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, N))
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("ww%02.0f-L%02.0f-N%i-hc.pdf" % (ww*100, L, N))

plt.figure('s')
plt.title('Configurational entropy for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, N))
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("ww%02.0f-L%02.0f-N%i-S.pdf" % (ww*100, L, N))