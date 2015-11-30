#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles
import readandcompute

if len(sys.argv) != 5:
    print 'useage: %s i ww L N ' % sys.argv[0]
    exit(1)

R = 1 # sphere radius

i = float(sys.argv[1])
#arg i = [0,1,2]

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[3])
#arg L = [5]

Ns = eval(sys.argv[4])
#arg Ns = [range(0,10)]
#arg N = 10

#seed = int(sys.argv[5])
#arg seed = [0]


Etas = map(lambda N: N*(4*np.pi/3*R**3/(L*2**i)**3), Ns) # works for i = 0,1

for N in Ns:
    T, U, F, CV, S, min_T = readandcompute.T_u_F_cv_s_minT('../data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, N))
    eta = Etas[Ns.index(N)]
    # plot multiple F(T) at fixed i
    plt.figure('F-Fix_i')
    plt.plot(F[750], eta, label="$T=%04.1f$" % T[750])
    plt.plot(F[250], eta, label="$T=%04.1f$" % T[250])
    
plt.figure('F-Fix_i')
plt.title('Absolute free energies for $\lambda=%g$, $L=%g$ and $i=%d$' % (ww, L, i))
plt.xlabel('$\eta$')
plt.ylabel('$F/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("ww%02.0f-L%02.0f-F_i%01d_fixed_i.pdf" % (ww*100, L, i))
