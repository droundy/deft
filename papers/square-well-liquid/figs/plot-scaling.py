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
#arg versions = [["nw", "wang_landau", "gaussian", "flat", "walkers", "kT2", "kT1"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-rt.dat" % (ww, ff, N, version) for version in versions for N in [5,6,7]]

# input: "data/periodic-ww*-ff*-N*-*-rt.dat"

N_regex = re.compile(r'-N([0-9]+)')
initialization_iters_regex = re.compile(r'# iterations:\s+([0-9]+)')

plt.title('Scaling for $\lambda=%g$, $\eta=%g$' % (ww, ff))
for version in versions:
    for filename in glob.glob("data/periodic-ww%04.2f-ff%04.2f-N*-%s-rt.dat" % (ww, ff, version)):
        data = numpy.loadtxt(filename)

        N = int(N_regex.findall(filename)[0])

        with open(filename, 'r') as content_file:
            content = content_file.read()

        init_iters = int(initialization_iters_regex.findall(content)[0])
        print 'N is', N, 'version is', version, 'iters', init_iters

        if N == 6:
            plt.semilogy(N, init_iters, styles.color[version]+'.', label=version)
        else:
            plt.semilogy(N, init_iters, styles.color[version]+'.')

        e_hist = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version))
        lnw = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version))
        energy = -e_hist[:,0]/N
        log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
        log10_dos = numpy.log10(e_hist[:,1]) - log10w
        log10_dos -= log10_dos.max()

plt.xlabel('$N$')
plt.ylabel('initialization iterations')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling.pdf" % (ww*100, ff*100))

if 'show' in sys.argv:
    plt.show()
