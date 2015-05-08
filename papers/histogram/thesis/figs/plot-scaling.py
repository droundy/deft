#!/usr/bin/env python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, re, string, os

sys.path.append('../figs')
import styles
import readandcompute

thesis_dir = os.getcwd()
os.chdir('../')

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

if len(sys.argv) != 7:
    print 'useage: %s ww ff Ns methods golden seed' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

Ns = eval(sys.argv[3])
#arg Ns = [range(5,26)]

methods = eval(sys.argv[4])
#arg methods = [["cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

golden = sys.argv[5]
#arg golden = ['tmmc-golden']

seed = int(sys.argv[6])
#arg seed = [0]

u_errors = {}
cv_errors = {}
s_errors = {}

for N in Ns:
    u_cv_s = readandcompute.u_cv_s(ww, ff, N, golden)
    if u_cv_s != None:
        for method in methods:
            u_cv_s_method = readandcompute.u_cv_s(ww, ff, N, method, seed)
            if u_cv_s_method != None:
                du = abs(u_cv_s_method[0]-u_cv_s[0]).max()
                if method not in u_errors:
                    u_errors[method] = numpy.array([[N, du]])
                else:
                    u_errors[method] = numpy.vstack([u_errors[method], [N, du]])
                dcv = abs(u_cv_s_method[1]-u_cv_s[1]).max()
                if method not in cv_errors:
                    cv_errors[method] = numpy.array([[N, dcv]])
                else:
                    cv_errors[method] = numpy.vstack([u_errors[method], [N, dcv]])
                ds = abs(u_cv_s_method[2]-u_cv_s[2]).max()
                if method not in s_errors:
                    s_errors[method] = numpy.array([[N, ds]])
                else:
                    s_errors[method] = numpy.vstack([u_errors[method], [N, ds]])

os.chdir(thesis_dir)

min_T = 0.2 # FIXME maybe shouldn't hardcode this?

plt.figure(figsize=(6,5))
for method in methods:
    plt.semilogy(u_errors[method][:,0], u_errors[method][:,1],
                 '.'+styles.plot(method),label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/u-scaling.pdf")

plt.figure(figsize=(6,5))
for method in methods:
    plt.semilogy(cv_errors[method][:,0], cv_errors[method][:,1],
                 '.'+styles.plot(method),label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel('$\\Delta C_V/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/cv-scaling.pdf")

plt.figure(figsize=(6,5))
for method in methods:
    plt.semilogy(s_errors[method][:,0], s_errors[method][:,1],
                 '.'+styles.plot(method),label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/s-scaling.pdf")
