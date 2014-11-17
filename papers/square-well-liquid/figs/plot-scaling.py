#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, glob, re
import styles

if len(sys.argv) not in [4,5]:
    print 'useage: %s ww ff N versions show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

versions = eval(sys.argv[3])
#arg versions = [["wang_landau", "gaussian", "flat", "walkers"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-rt.dat" % (ww, ff, N, version) for version in versions for N in [5,6,7,8,9]]

# input: "data/periodic-ww*-ff*-N*-*-rt.dat"

N_regex = re.compile(r'-N([0-9]+)')
initialization_iters_regex = re.compile(r'# iterations:\s+([0-9]+)')

plt.title('Scaling for $\lambda=%g$, $\eta=%g$' % (ww, ff))
Ns = {}
init_iters = {}
for version in versions:
    Ns[version] = []
    init_iters[version] = []
    for filename in glob.glob("data/periodic-ww%04.2f-ff%04.2f-N*-%s-rt.dat" % (ww, ff, version)):
        data = numpy.loadtxt(filename)

        N = int(N_regex.findall(filename)[0])
        Ns[version].append(N)

        with open(filename, 'r') as content_file:
            content = content_file.read()

        init_iters[version].append(int(initialization_iters_regex.findall(content)[0]))

        e_hist = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version))
        lnw = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version))
        energy = -e_hist[:,0]/N
        log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
        log10_dos = numpy.log10(e_hist[:,1]) - log10w
        log10_dos -= log10_dos.max()

    # The following sorts the lists according to the number of spheres
    indexes = range(len(Ns[version]))
    indexes.sort(key=Ns[version].__getitem__)
    Ns[version] = map(Ns[version].__getitem__, indexes)
    init_iters[version] = map(init_iters[version].__getitem__, indexes)

    plt.semilogy(Ns[version], init_iters[version], styles.color[version]+'.-', label=version)

plt.xlabel('$N$')
plt.ylabel('initialization iterations')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling.pdf" % (ww*100, ff*100))

if 'show' in sys.argv:
    plt.show()
