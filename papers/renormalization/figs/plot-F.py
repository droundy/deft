#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# This code potentially needs 

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles # nice plots
import readandcompute # needed for many functions, look at this if you haven't!

if len(sys.argv) != 5:
    print 'useage: %s i ww L [Ns] RECURSION-LEVEL WELL-WIDTH LENGTH [SPHERES]' % sys.argv[0]
    exit(1)

# Necessary constants and parameters 

R = 1 # sphere radius

i = float(sys.argv[1])
#arg i = [0,1,2]

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(sys.argv[3])
#arg L = [5]

#seed = int(sys.argv[5])
#arg seed = [0]

Ts = [0.2,0.5,0.8,1.5,2.0]
max_T = Ts[-1]
# Fixed T for plots


styles = ['ro','bo','g-','c-']

plot_all = plot_F_T = plot_F_eta = False

if len(sys.argv) < 5:
    # If no Ns are given, plot everything we have
    plot_all = True
else:
    # Choose plot type here. Might make this an argument, but seems unnecessary
    plot_F_T = True
    Ns = eval(sys.argv[4])
    #arg Ns = [range(0,10)]
 
    
if plot_all:
    # Is this necessary?
    for l in range(0, len(Ts)):
        eta, U, F, CV, S, min_T = readandcompute.eta_u_F_cv_s_minT('data/scrunched-ww%04.2f-L%04.2f' % (ww, L), T, max_T)
        plt.figure('F-all')
        plt.plot(eta, F, styles[j], label="$T=%02.1f$" % T)
        plt.xlabel('$\eta$')
        output_filename = "figs/ww%02.0f-L%02.0f-F_i%01d_F-eta_all.pdf" % (ww*100, L, i)
elif plot_F_eta:
    for N in set(Ns):
        eta = N*(4*np.pi/3*R**3/(L*2**i)**3)
        for j in range(0,len(Ts)):            
            T, U, F, CV, S, min_T = readandcompute.T_u_F_cv_s_minT('data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, N), max_T)
            q = readandcompute.nearest_T(T,Ts[j])
            Ts[j] = T[q]

            plt.figure('F-Fix_i')
            if N == Ns[0]:
                #hacky label magic
                plt.plot(eta, F[q], styles[j%len(styles)], label="$T=%02.1f$" % Ts[j])
            else:
                plt.plot(eta, F[q], styles[j%len(styles)])
                
    output_filename = "figs/ww%02.0f-L%02.0f-i%01d_F-eta.pdf" % (ww*100, L, i)
    plt.figure('F-eta')
    plt.xlabel('$\eta$')
    
elif plot_F_T:
    etas = np.sort(map(lambda N: N*(4*np.pi/3*R**3/(L*2**i)**3), Ns))
    for k in range(0,len(etas)):
        T, U, F, CV, S, min_T = readandcompute.T_u_F_cv_s_minT('data/scrunched-ww%04.2f-L%04.2f/i%01d/N%03d/data' % (ww, L,  i, Ns[k]), max_T)
        plt.figure('F-T')
        plt.plot(T, F, styles[k%len(styles)],label= "$\eta=%01.2f$" % etas[k])
    plt.xlabel('$T/\epsilon$')
    plt.axvline(x=min_T)
    output_filename = "figs/ww%02.0f-L%02.0f-i%01d_F-T.pdf" % (ww*100, L, i)
        
plt.title('Absolute free energies for $\lambda=%g$, $L=%g$ and $i=%d.$' % (ww, L, i))

plt.ylabel('$F/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)


plt.savefig(output_filename)



