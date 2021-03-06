#!/usr/bin/env python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, re, string, os
import styles
import readandcompute

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

if len(sys.argv) != 6:
    print('useage: %s ww ff methods golden seed' % sys.argv[0])
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

methods = eval(sys.argv[3])
#arg methods = [["simple_flat","tmmc","oetmmc"]]

golden = sys.argv[4]
#arg golden = ['tmmc-golden']

seed = int(sys.argv[5])
#arg seed = [0]

all_Ns = os.popen("ls data/periodic-ww%.2f-ff%.2f.*-golden-E.dat"
                  %(ww, ff)).readlines()
all_Ns = numpy.sort([ int(N.split('-')[-4][1:]) for N in all_Ns ])

if len(all_Ns) == 0:
    all_Ns = [5, 10, 20]

######################################################################
# initialization scaling info
######################################################################

initialization_iters_regex = re.compile(r'# iterations:\s+([0-9]+)')

plt.title('Scaling for $\lambda=%g$, $\eta=%g$' % (ww, ff))
init_iters = {}
Emins = {}
samples = {}
for method in methods:
  init_iters[method] = []
  Emins[method] = []
  samples[method] = []

  for N in all_Ns:
    filename = "data/s%03d/periodic-ww%04.2f-ff%04.2f-N%d-%s-lnw.dat" \
               % (seed, ww, ff, N, method)
    wildfilename = "data/s%03d/periodic-ww%04.2f-ff%04.2f-N%d-%s-%%s.dat" \
                   % (seed, ww, ff, N, method)

    with open(filename, 'r') as content_file:
        content = content_file.read()
    init_iters[method].append(int(initialization_iters_regex.findall(content)[0]))

    E_data = numpy.loadtxt(wildfilename % 'E', ndmin=2)
    Emins[method].append(E_data[:, 0].max())

    sample_data = numpy.loadtxt(wildfilename % 'ps', ndmin=2)
    samples[method].append(sample_data[len(sample_data[:, 1])-1, 1])

  plt.figure('iters')
  plt.semilogy(all_Ns, init_iters[method], styles.color(method)+'.-',
               label=styles.title(method))
  plt.figure('emin')
  plt.plot(all_Ns, Emins[method], styles.color(method)+'.-', label=styles.title(method))

plt.figure('iters')
plt.xlabel('$N$')
plt.ylabel('Initialization iterations')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling.pdf" % (ww*100, ff*100))

plt.figure('emin')
plt.xlabel('$N$')
plt.ylabel('Emin')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling-emin.pdf" % (ww*100, ff*100))

######################################################################
# scaling info table
######################################################################

texmethods = [m.replace('_', ' ') for m in methods]

tex = open("figs/scaling-table-ww%02.0f-ff%02.0f.tex" % (ww*100, ff*100), "w")
tex.write(r"""
\begin{tabular}{c|%s}
""" % ("c"*len(methods)))
for method in texmethods:
    tex.write(" & " +method)
tex.write(r"""\\
\hline\hline
""")

for i in range(len(all_Ns)):
    N = all_Ns[i]
    tex.write(r""" N = %d \\
  initialization""" % N)

    tex.write(string.join([" & %d " % init_iters[m][i] for m in methods]))
    tex.write("\\\\\n")

    tex.write('Emin')
    tex.write(string.join([" & %d " % Emins[m][i] for m in methods]))
    tex.write("\\\\\n")

    tex.write('samples')
    tex.write(string.join([" & %d " % samples[m][i] for m in methods]))
    tex.write("\\\\\n")

tex.write(r"""\end{tabular}
""")

######################################################################
# error scaling figures
######################################################################

u_errors = {}
cv_errors = {}
s_errors = {}

for N in all_Ns:
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

min_T = 0.2 # FIXME maybe shouldn't hardcode this?

plt.figure()
for method in list(u_errors.keys()):
    plt.semilogy(u_errors[method][:, 0], u_errors[method][:, 1],
                 '.'+styles.plot(method), label=styles.title(method))

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel('$N$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-u_errors.pdf" % (ww*100, ff*100))

plt.figure()
for method in list(cv_errors.keys()):
    plt.semilogy(cv_errors[method][:, 0], cv_errors[method][:, 1],
                 '.'+styles.plot(method), label=styles.title(method))

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel('$N$')
plt.ylabel('$\\Delta C_V/N\epsilon$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-cv_errors.pdf" % (ww*100, ff*100))

plt.figure()
for method in list(s_errors.keys()):
    plt.semilogy(s_errors[method][:, 0], s_errors[method][:, 1],
                 '.'+styles.plot(method), label=styles.title(method))

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel('$N$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-s_errors.pdf" % (ww*100, ff*100))

######################################################################
# error vs reference error figures
######################################################################

if 'cfw' in methods:
    reference = 'cfw'
else:
    reference = 'tmmc'

comp_methods = list(s_errors.keys())
comp_methods.remove(reference)

plt.figure()
for method in comp_methods:
    plt.loglog(u_errors[reference][:, 1], u_errors[method][:, 1],
               '.'+styles.color(method), label=styles.title(method))
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.plot([xmin, xmax], [ymin, ymax], 'k:')

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel('%s $\\Delta U/N\epsilon$'%styles.title(reference))
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-u_error_comp.pdf" % (ww*100, ff*100))

plt.figure()
for method in comp_methods:
    plt.loglog(cv_errors[reference][:, 1], cv_errors[method][:, 1],
               '.'+styles.color(method), label=styles.title(method))
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.plot([xmin, xmax], [ymin, ymax], 'k:')

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel('%s $\\Delta C_V/N\epsilon$'%styles.title(reference))
plt.ylabel('$\\Delta C_V/N\epsilon$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-cv_error_comp.pdf" % (ww*100, ff*100))

plt.figure()
for method in comp_methods:
    plt.loglog(s_errors[reference][:, 1], s_errors[method][:, 1],
               '.'+styles.color(method), label=styles.title(method))
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.plot([xmin, xmax], [ymin, ymax], 'k:')

plt.title('Maximum error with $\lambda=%g$, $\eta=%g$, and $T_{min}=%g$' % (ww, ff, min_T))
plt.xlabel(r'%s $\Delta S_{\textit{config}}/Nk$'%styles.title(reference))
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-s_error_comp.pdf" % (ww*100, ff*100))
