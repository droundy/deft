#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, glob, re, string
import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N methods show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

all_Ns = eval(sys.argv[3])
#arg all_Ns = [[5,7,10,20]]

methods = eval(sys.argv[4])
#arg methods = [["wang_landau","simple_flat","optimized_ensemble","tmmc", "oetmmc"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, method, dat) for method in methods for N in all_Ns for dat in ['ps', 'lnw', 'E']]

N_regex = re.compile(r'-N([0-9]+)')
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
    filename = "data/periodic-ww%04.2f-ff%04.2f-N%d-%s-lnw.dat" % (ww, ff, N, method)
    wildfilename = "data/periodic-ww%04.2f-ff%04.2f-N%d-%s-%%s.dat" % (ww, ff, N, method)

    with open(filename, 'r') as content_file:
        content = content_file.read()
    init_iters[method].append(int(initialization_iters_regex.findall(content)[0]))

    E_data = numpy.loadtxt(wildfilename % 'E', ndmin=2)
    Emins[method].append(E_data[:, 0].max())

    sample_data = numpy.loadtxt(wildfilename % 'ps', ndmin=2)
    samples[method].append(sample_data[len(sample_data[:,1])-1, 1])

  plt.figure(1)
  plt.semilogy(all_Ns, init_iters[method], styles.color(method)+'.-', label=method)
  plt.figure(2)
  plt.plot(all_Ns, Emins[method], styles.color(method)+'.-', label=method)

plt.figure(1)
plt.xlabel('$N$')
plt.ylabel('Initialization iterations')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling.pdf" % (ww*100, ff*100))

plt.figure(2)
plt.xlabel('$N$')
plt.ylabel('Emin')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling-emin.pdf" % (ww*100, ff*100))

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

