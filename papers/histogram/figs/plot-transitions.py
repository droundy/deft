#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import styles

if len(sys.argv) != 7:
    print 'useage: %s ww ff N seed' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = range(5,21)+[100, 200, 1000]

golden_method = sys.argv[4]
#arg golden_method = ['tmmc-golden']

method = sys.argv[5]
#arg method = ['tmmc']

seed = int(sys.argv[6])
#arg seed = [0]

bothdata = [0,0]

golden_data_file = "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-transitions.dat" \
                   % (ww, ff, N, golden_method)
method_data_file = "data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-transitions.dat" \
                   % (seed, ww, ff, N, method)

golden_data = numpy.loadtxt(golden_data_file)
method_data = numpy.loadtxt(method_data_file)

def get_de(data_file):
    with open(data_file) as f:
        for line in f:
            if '# energy\t' in line:
                de = [ -int(val) for val in line.split()[2:] ]
                break
    return de

de_golden = get_de(golden_data_file)
de_method = get_de(method_data_file)

for method, data, de in [(golden_method, golden_data, de_golden), (method, method_data, de_method)]:
    plt.figure()
    plt.title('Energy transitions from %s for $\lambda=%g$, $n^*=%g$, and $N=%i$' % (method, ww, ff, N))

    energy = -data[:,0]
    maxde = 0
    minde = 0
    transitions = data[:,1:]
    de, e = numpy.meshgrid(de, energy)

    for i in range(len(energy)):
        norm = transitions[i,:].sum()
        if norm > 0:
            transitions[i,:] /= norm
        for j in range(len(transitions[i,:])):
            if transitions[i,j] > 0:
                if de[i,j] < minde:
                    minde = de[i,j]
                if de[i,j] > maxde:
                    maxde = de[i,j]
    c = plt.pcolor(de+0.5, e, transitions, vmin=0, vmax=1)
    plt.colorbar(c)
    plt.xlim(minde,maxde)
    plt.ylim(e.min(), e.max())

    plt.xlabel('$\Delta E/\epsilon$')
    plt.ylabel('$E/\epsilon$')

    if method == golden_method:
        plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-transitions.pdf" % (ww*100, ff*100, N, golden_method))
    if method == method:
        plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-transitions.pdf" % (ww*100, ff*100, N, method))

plt.figure()

if golden_data[0,0] < method_data[0,0]:
    golden_data = golden_data[method_data[0,0] - golden_data[0,0]:,:]
else:
    method_data = method_data[golden_data[0,0] - method_data[0,0]:,:]

energy = -golden_data[:,0]

# fixme: do reasonable things when de_golden != de_method ?
de, e = numpy.meshgrid(de_golden, energy)

cut2 = golden_data[:,0].min() - method_data[:,0].min()

for i in range(len(golden_data[:,0])):
    if i < len(method_data[:,0]):
        for j in range(len(golden_data[1,:])):
            golden_data[i,j] -= method_data[i,j]

c = plt.pcolor(de-0.5, e, golden_data[:,1:], cmap='RdBu', vmin=-0.05, vmax = 0.05)
plt.colorbar(c)
plt.xlim(-8,8)
plt.ylim(e.min(), e.max())
plt.xlabel('$\Delta E/\epsilon$')
plt.ylabel('$E/\epsilon$')

plt.title('TM discrepancy between %s and %s\nfor $\lambda=%g$, $n^*=%g$, and $N=%i$' % (golden_method, method, ww, ff, N))

plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-%s-compare-transitions.pdf" % (ww*100, ff*100, N, golden_method, method))

