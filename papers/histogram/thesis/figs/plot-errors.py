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

if len(sys.argv) != 8:
    print 'useage: %s ww ff Ns seeds methods golden reference' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

Ns = eval(sys.argv[3])
#arg Ns = [range(5,26)]

seeds = eval(sys.argv[4])
#arg seeds = [range(30)]

methods = eval(sys.argv[5])
#arg methods = [["simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

golden = sys.argv[6]
#arg golden = ['tmmc-golden']

reference = sys.argv[7]
#arg reference = ['cfw']

fig_u = plt.figure('u',figsize=(6.3,4))
ax_u = fig_u.add_axes([0.12, 0.12, 0.7, 0.85])

fig_cv = plt.figure('cv',figsize=(6.3,4))
ax_cv = fig_cv.add_axes([0.1, 0.13, 0.7, 0.83])

fig_s = plt.figure('s',figsize=(6.3,4))
ax_s = fig_s.add_axes([0.1, 0.13, 0.73, 0.85])

figs = ['u','cv','s']

Tmax = 0.2

def mean_errs(ww, ff, N, method, golden_data):
    u_err = 0.0
    cv_err = 0.0
    s_err = 0.0
    numseeds = 0
    for seed in seeds:
        t_u_cv_s_method = readandcompute.t_u_cv_s(ww, ff, N, method, seed)
        if t_u_cv_s_method != None:
            du = abs((t_u_cv_s_method[1]-golden_data[1])[t_u_cv_s_method[0] > Tmax]).max()
            dcv = abs((t_u_cv_s_method[2]-golden_data[2])[t_u_cv_s_method[0] > Tmax]).max()
            ds = abs((t_u_cv_s_method[3]-golden_data[3])[t_u_cv_s_method[0] > Tmax]).max()
            u_err += du
            cv_err += dcv
            s_err += ds
            numseeds += 1
    if numseeds > 0:
        return numpy.array([u_err,cv_err,s_err])/numseeds/N
    return None

for N in Ns:
    t_u_cv_s = readandcompute.t_u_cv_s(ww, ff, N, golden)
    golden_data = t_u_cv_s
    reference_error = mean_errs(ww, ff, N, reference, golden_data)
    print ' N',N
    for method in methods:
        method_err = mean_errs(ww, ff, N, method, golden_data)
        print '  ',method
        print '   ',method_err
        for i in range(len(figs)):
            plt.figure(figs[i])
            plt.loglog(reference_error[i], method_err[i], styles.dots(method),
                       markerfacecolor='none', markeredgecolor=styles.color(method),
                       label=styles.title(method))

os.chdir(thesis_dir)

for fig in figs:
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
plt.legend(handle_list, label_list,loc='lower right',bbox_to_anchor=(1.25,0))

plt.xlabel(r'%s $\Delta U/N\epsilon$'%styles.title(reference))
plt.ylabel(r'$\Delta U/N\epsilon$')
plt.savefig("figs/internal-energy-comps.pdf")

plt.figure('cv')
plt.legend(handle_list, label_list,loc='lower right',bbox_to_anchor=(1.25,0))
plt.xlabel(r'%s $\Delta C_V/Nk$'%styles.title(reference))
plt.ylabel(r'$\Delta C_V/Nk$')
plt.savefig("figs/heat-capacity-comps.pdf")

plt.figure('s')
plt.legend(handle_list, label_list,loc='lower right',bbox_to_anchor=(1.25,0))
plt.xlabel(r'%s $\Delta S_{\textrm{config}}/Nk$'%styles.title(reference))
plt.ylabel(r'$\Delta S_{\textrm{config}}/Nk$')
plt.savefig("figs/config-entropy-comps.pdf")

