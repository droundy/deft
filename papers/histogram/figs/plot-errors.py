#!/usr/bin/env python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, re, string, os
import styles
import readandcompute

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

all_goldens = os.popen("git ls-files | grep 'data/periodic-ww.*-golden-E.dat'").readlines()

golden = 'tmmc-golden'

methods = ['tmmc', 'oetmmc', 'simple_flat', 'wang_landau', 'vanilla_wang_landau']
reference = 'cfw'
seeds = range(0, 30)

######################################################################
# Make beautiful error plot
######################################################################

plt.figure(1)

Tmax = 0.2

def mean_u_err(ww, ff, N, method, golden_u):
    toterr = 0.0
    numseeds = 0
    for seed in seeds:
        t_u_cv_s_method = readandcompute.t_u_cv_s(ww, ff, N, method, seed)
        if t_u_cv_s_method != None:
            du = abs((t_u_cv_s_method[1]-golden_u)[t_u_cv_s_method[0] > Tmax]).max()
            lw = 0.01
            if seed == 0:
                plt.plot(t_u_cv_s_method[0], (t_u_cv_s_method[1]-golden_u)/N,
                         styles.plot(method), label=styles.title(method), linewidth=lw)
            else:
                plt.plot(t_u_cv_s_method[0], (t_u_cv_s_method[1]-golden_u)/N,
                         styles.plot(method), linewidth=lw)
            toterr += du
            numseeds += 1
    if numseeds > 0:
        return toterr/numseeds/N
    return None

labels_added = set()

for golden_comp in all_goldens:
    bits = golden_comp.split('-')
    ww = float(bits[1][2:])
    ff = float(bits[2][2:])
    N = int(bits[3][1:])
    print 'ww = %g, ff = %g, N = %d' % (ww, ff, N)
    t_u_cv_s = readandcompute.t_u_cv_s(ww, ff, N, golden)
    if t_u_cv_s != None:
        plt.figure(2)
        plt.cla()
        plt.title(r'error in $U/N$ for $N=%d$ and $\eta=%g$ and $\lambda=%g$' % (N, ff, ww))
        plt.xlabel(r'$T$')
        plt.ylabel(r'$\Delta U$')

        golden_u = t_u_cv_s[1]
        reference_error = mean_u_err(ww, ff, N, reference, golden_u)
        if reference_error == None:
            continue
        print '    ', reference, reference_error
        for method in methods:
            plt.figure(2)
            method_err = mean_u_err(ww, ff, N, method, golden_u)
            if method_err != None:
                print '    ', method, method_err
                plt.figure(1)
                if method not in labels_added:
                    plt.loglog(reference_error, method_err, styles.dots(method),
                               markerfacecolor='none', markeredgecolor=styles.color(method),
                               label=styles.title(method))
                    labels_added |= set([method])
                else:
                    plt.loglog(reference_error, method_err, styles.dots(method),
                               markerfacecolor='none', markeredgecolor=styles.color(method))
                plt.figure(2)
        plt.figure(2)
        plt.legend(loc='best')
        plt.savefig(
            'figs/energy-for-debugging-ww%04.2f-ff%04.2f-N%i.pdf' % (ww, ff, N))

min_T = 0.2 # FIXME maybe shouldn't hardcode this?

plt.figure(1)
xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.axis('equal')
plt.loglog([min(xmin,ymin),max(xmax,ymax)],[min(xmin,ymin),max(xmax,ymax)],'k:')
plt.legend(loc='best')

plt.title('average maximum error in $U$ for various configurations')
plt.xlabel(r'%s $\left<\Delta U\right>/N\epsilon$'%styles.title(reference))
plt.ylabel(r'$\Delta U/N\epsilon$')
#plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/u-error-comparison.pdf")
