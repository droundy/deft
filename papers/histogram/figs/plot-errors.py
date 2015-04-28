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

plt.figure()

def mean_u_err(ww, ff, N, method, golden_u):
    toterr = 0.0
    numseeds = 0
    for seed in seeds:
        u_cv_s_method = readandcompute.u_cv_s(ww, ff, N, method, seed)
        if u_cv_s_method != None:
            du = abs(u_cv_s_method[0]-golden_u).max()
            toterr += du
            numseeds += 1
    if numseeds > 0:
        return toterr/numseeds/N
    return None

for golden_comp in all_goldens:
    bits = golden_comp.split('-')
    ww = float(bits[1][2:])
    ff = float(bits[2][2:])
    N = int(bits[3][1:])
    print 'ww = %g, ff = %g, N = %d' % (ww, ff, N)
    u_cv_s = readandcompute.u_cv_s(ww, ff, N, golden)
    if u_cv_s != None:
        golden_u = u_cv_s[0]
        reference_error = mean_u_err(ww, ff, N, reference, golden_u)
        if reference_error == None:
            continue
        print '    ', reference, reference_error
        for method in methods:
            method_err = mean_u_err(ww, ff, N, method, golden_u)
            if method_err != None:
                print '    ', method, method_err
                plt.plot(reference_error, method_err, styles.dots(method))

min_T = 0.2 # FIXME maybe shouldn't hardcode this?

xmin,xmax = plt.xlim()
ymin,ymax = plt.ylim()
plt.plot([min(xmin,ymin),max(xmax,ymax)],[min(xmin,ymin),max(xmax,ymax)],'k:')

plt.title('average maximum error in $U$ for various configurations')
plt.xlabel('%s $\\Delta U/N\epsilon$'%styles.title(reference))
plt.ylabel('$\\Delta U/N\epsilon$')
#plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/u-error-comparison.pdf")
