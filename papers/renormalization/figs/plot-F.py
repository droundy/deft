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
    print 'useage: %s i ww L [Ns] RECURSION-LEVEL WELL-WIDTH LENGTH [SPHERES]' % sys.argv[0]
    exit(1)

R = 1 # sphere radius

i = float(sys.argv[1])
#arg i = [0,1,2]

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[3])
#arg L = [5]

#seed = int(sys.argv[5])
#arg seed = [0]

Ts = [0.1,0.3,0.5,0.8,1.5]
# Fixed T for plots

#find closest Ts in T_range


styles = ['ro','bo','rx','bx','gx']

plot_all = False

if len(sys.argv) < 5:
    plot_all = True
else:
    Ns = eval(sys.argv[4])
    #arg Ns = [range(0,10)]

    
if plot_all:
    for T in set(Ts):
        eta, U, F, CV, S, min_T = readandcompute.eta_u_F_cv_s_minT('data/scrunched-ww%04.2f-L%04.2f' % (ww, L), T)
        plt.figure('F-all')
        plt.plot(eta, F, styles[j], label="$T=%02.1f$" % T)
        
else:
    for N in set(Ns):
        eta = N*(4*np.pi/3*R**3/(L*2**i)**3)
        for j in range(0,len(Ts)):            
            T, U, F, CV, S, min_T = readandcompute.T_u_F_cv_s_minT('data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, N))
            q = readandcompute.nearest_T(T,Ts[j])
            Ts[j] = T[q]

            plt.figure('F-Fix_i')
            if N == Ns[0]:
                #hacky label magic
                plt.plot(eta, F[q], styles[j%5], label="$T=%02.1f$" % Ts[j])
            else:
                plt.plot(eta, F[q], styles[j%5])
    plt.figure('F-Fix_i')
                 
plt.title('Absolute free energies for $\lambda=%g$, $L=%g$ and $i=%d.$' % (ww, L, i))
plt.xlabel('$\eta$')
plt.ylabel('$F/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
output_filename = "figs/ww%02.0f-L%02.0f-F_i%01d_fixed_i.pdf" % (ww*100, L, i)
plt.savefig(output_filename)



