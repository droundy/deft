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
#arg i = [0,1,2]

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[3])
#arg L = [5]

Ns = eval(sys.argv[4])
#arg Ns = [range(0,10)]
#arg N = 10

method = "tmmc"

#seed = int(sys.argv[5])
#arg seed = [0]

# for method in set(methods):
#     T, U, F, CV, S, min_T = readandcompute.T_u_f_cv_s_minT('../data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, N))

    
#     S -= readandcompute.absolute_f('../data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/absolute/' % (ww, L,  i, N))
    
#     plt.figure('u')
#     plt.plot(T,U/N,styles.plot(method),label=styles.title(method))

#     plt.figure('hc')
#     plt.plot(T,CV/N,styles.plot(method),label=styles.title(method))

#     plt.figure('s')
#     plt.plot(T,S/N,styles.plot(method),label=styles.title(method))
for N in Ns:
        T, U, F, CV, S, min_T = readandcompute.T_u_F_cv_s_minT('../data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, N))
        plt.figure('F-T')
        plt.plot(F, T,  label="$N=%d$" % N)
    
plt.figure('F-T')
plt.title('Absolute free energies for $\lambda=%g$, $L=%g$ and $i=%d$' % (ww, L, i))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$F/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("ww%02.0f-L%02.0f-i%i-F_eta.pdf" % (ww*100, L, i))

plt.figure('F-eta')
plt.title('Absolute free energies for $\lambda=%g$, $L=%g$, and $i=%d$' % (ww,L,i))
plt.xlabel('$\eta$')
plt.ylabel('$F/N\epsilon$')
plt.legend(loc='best')


# plt.figure('u')
# plt.title('Specific internal energy for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, n))
# plt.xlabel('$kT/\epsilon$')
# plt.ylabel('$U/N\epsilon$')
# plt.legend(loc='best')
# plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
# plt.tight_layout(pad=0.2)
# plt.savefig("ww%02.0f-L%02.0f-N%i-u.pdf" % (ww*100, L, N))

# plt.figure('hc')
# plt.title('Specific heat capacity for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, N))
# plt.ylim(0)
# plt.xlabel('$kT/\epsilon$')
# plt.ylabel('$C_V/Nk$')
# plt.legend(loc='best')
# plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
# plt.tight_layout(pad=0.2)
# plt.savefig("ww%02.0f-L%02.0f-N%i-hc.pdf" % (ww*100, L, N))

# plt.figure('s')
# plt.title('Configurational entropy for $\lambda=%g$, $L=%g$, and $N=%i$' % (ww, L, N))
# plt.xlabel(r'$kT/\epsilon$')
# plt.ylabel(r'$S_{\textit{config}}/Nk$')
# plt.legend(loc='best')
# plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
# plt.tight_layout(pad=0.2)
# plt.savefig("ww%02.0f-L%02.0f-N%i-S.pdf" % (ww*100, L, N))
