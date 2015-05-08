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

if len(sys.argv) != 9:
    print 'useage: %s ww ff Ns min_T methods golden reference seeds' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

Ns = eval(sys.argv[3])
#arg Ns = [range(5,26)]

min_T = float(sys.argv[4])
#arg min_T = [0.2]

methods = eval(sys.argv[5])
#arg methods = [["cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

golden = sys.argv[6]
#arg golden = ['tmmc-golden']

reference = sys.argv[7]
#arg reference = ['cfw']

seeds = eval(sys.argv[8])
#arg seeds = [range(10,21)]

u_errors = {}
cv_errors = {}
s_errors = {}
for method in methods:
    u_errors[method] = numpy.zeros(len(Ns))
    cv_errors[method] = numpy.zeros(len(Ns))
    s_errors[method] = numpy.zeros(len(Ns))

for n in range(len(Ns)):
    golden_data = readandcompute.t_u_cv_s(ww, ff, Ns[n], golden)
    for method in methods:
        for seed in seeds:
            t_u_cv_s_method = readandcompute.t_u_cv_s(ww, ff, Ns[n], method, seed)
            du = abs((t_u_cv_s_method[1]-golden_data[1])[t_u_cv_s_method[0] > min_T]).max()
            dcv = abs((t_u_cv_s_method[2]-golden_data[2])[t_u_cv_s_method[0] > min_T]).max()
            ds = abs((t_u_cv_s_method[3]-golden_data[3])[t_u_cv_s_method[0] > min_T]).max()

            u_errors[method][n] += du
            cv_errors[method][n] += dcv
            s_errors[method][n] += ds

        u_errors[method][n] /= len(seeds)
        cv_errors[method][n] /= len(seeds)
        s_errors[method][n] /= len(seeds)

os.chdir(thesis_dir)


### scaling figures
fig_size = (6,5)

plt.figure(figsize=fig_size)
for method in methods:
    plt.semilogy(Ns, u_errors[method],'.'+styles.plot(method),
                 label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/u-scaling.pdf")


plt.figure(figsize=fig_size)
for method in methods:
    plt.semilogy(Ns, cv_errors[method],'.'+styles.plot(method),
                 label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel('$\\Delta C_V/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/cv-scaling.pdf")

plt.figure(figsize=fig_size)
for method in methods:
    plt.semilogy(Ns, s_errors[method],'.'+styles.plot(method),
                 label=styles.title(method))

plt.xlabel('$N$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/s-scaling.pdf")

### error vs error figures

methods.remove(reference)

plt.figure('u',figsize=fig_size)
plt.figure('cv',figsize=fig_size)
plt.figure('s',figsize=fig_size)

for n in range(len(Ns)):
    for method in methods:
        plt.figure('u')
        plt.loglog(u_errors[reference][n], u_errors[method][n], styles.dots(method),
                   markerfacecolor='none', markeredgecolor=styles.color(method),
                   label=styles.title(method).replace('Vanilla Wang-Landau','Wang-Landau'))
        plt.figure('cv')
        plt.loglog(cv_errors[reference][n], cv_errors[method][n], styles.dots(method),
                   markerfacecolor='none', markeredgecolor=styles.color(method),
                   label=styles.title(method).replace('Vanilla Wang-Landau','Wang-Landau'))
        plt.figure('s')
        plt.loglog(s_errors[reference][n], s_errors[method][n], styles.dots(method),
                   markerfacecolor='none', markeredgecolor=styles.color(method),
                   label=styles.title(method).replace('Vanilla Wang-Landau','Wang-Landau'))

for fig in ['u','cv','s']:
    plt.figure(fig)
    xmin,xmax = plt.xlim()
    ymin,ymax = plt.ylim()
    plt.axis('equal')
    plt.loglog([min(xmin,ymin),max(xmax,ymax)],[min(xmin,ymin),max(xmax,ymax)],'k:')

plt.figure('u')

handles, labels = plt.gca().get_legend_handles_labels()
handle_list, label_list = [], []
for handle, label in zip(handles, labels):
    if label not in label_list:
        handle_list.append(handle)
        label_list.append(label)
plt.legend(handle_list, label_list,loc='lower right')

plt.xlabel(r'%s $\Delta U/N\epsilon$'%styles.title(reference))
plt.ylabel(r'$\Delta U/N\epsilon$')
plt.savefig("figs/u-comps.pdf")

plt.figure('cv')
plt.legend(handle_list, label_list,loc='lower right')
plt.xlabel(r'%s $\Delta C_V/Nk$'%styles.title(reference))
plt.ylabel(r'$\Delta C_V/Nk$')
plt.savefig("figs/cv-comps.pdf")

plt.figure('s')
plt.legend(handle_list, label_list,loc='lower right')
plt.xlabel(r'%s $\Delta S_{\textrm{config}}/Nk$'%styles.title(reference))
plt.ylabel(r'$\Delta S_{\textrm{config}}/Nk$')
plt.savefig("figs/s-comps.pdf")

