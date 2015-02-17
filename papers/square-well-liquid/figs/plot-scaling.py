#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, glob, re, string
import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N versions show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

all_Ns = eval(sys.argv[3])
#arg all_Ns = [[5,8,10,20]]

versions = eval(sys.argv[4])
#arg versions = [["wang_landau","robustly_optimistic","optimized_ensemble","tmmc","kT0.5","kT0.4"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, version, dat) for version in versions for N in all_Ns for dat in ['ps', 'lnw', 'E']]

N_regex = re.compile(r'-N([0-9]+)')
initialization_iters_regex = re.compile(r'# iterations:\s+([0-9]+)')

plt.title('Scaling for $\lambda=%g$, $\eta=%g$' % (ww, ff))
init_iters = {}
Emins = {}
samples = {}
for version in versions:
  init_iters[version] = []
  Emins[version] = []
  samples[version] = []
  for N in all_Ns:
    filename = "data/periodic-ww%04.2f-ff%04.2f-N%d-%s-lnw.dat" % (ww, ff, N, version)
    wildfilename = "data/periodic-ww%04.2f-ff%04.2f-N%d-%s-%%s.dat" % (ww, ff, N, version)

    with open(filename, 'r') as content_file:
        content = content_file.read()
    init_iters[version].append(int(initialization_iters_regex.findall(content)[0]))

    E_data = numpy.loadtxt(wildfilename % 'E', ndmin=2)
    Emins[version].append(E_data[:, 0].max())

    sample_data = numpy.loadtxt(wildfilename % 'ps', ndmin=2)
    samples[version].append(sample_data[len(sample_data[:,1])-1, 1])

  plt.figure(1)
  plt.semilogy(all_Ns, init_iters[version], styles.color(version)+'.-', label=version)
  plt.figure(2)
  plt.plot(all_Ns, Emins[version], styles.color(version)+'.-', label=version)

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

texversions = [v.replace('_', ' ') for v in versions]

tex = open("figs/scaling-table-ww%02.0f-ff%02.0f.tex" % (ww*100, ff*100), "w")
tex.write(r"""
\begin{tabular}{c|%s}
""" % ("c"*len(versions)))
for ver in texversions:
    tex.write(" & " +ver)
tex.write(r"""\\
\hline\hline
""")

for i in range(len(all_Ns)):
    N = all_Ns[i]
    tex.write(r""" N = %d \\
  initialization""" % N)

    tex.write(string.join([" & %d " % init_iters[v][i] for v in versions]))
    tex.write("\\\\\n")

    tex.write('Emin')
    tex.write(string.join([" & %d " % Emins[v][i] for v in versions]))
    tex.write("\\\\\n")

    tex.write('samples')
    tex.write(string.join([" & %d " % samples[v][i] for v in versions]))
    tex.write("\\\\\n")

tex.write(r"""\end{tabular}
""")

if 'show' in sys.argv:
    plt.show()
