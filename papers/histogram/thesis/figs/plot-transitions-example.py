#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig_size = (5,4)

if len(sys.argv) != 5:
    print 'useage: %s ww ff N method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = int(sys.argv[3])
#arg N = [25]

method = sys.argv[4]
#arg method = ['tmmc-golden']

data_file = "../data/periodic-ww%04.2f-ff%04.2f-N%i-%s-transitions.dat" %(ww, ff, N, method)

data = numpy.loadtxt(data_file)

def get_de(data_file):
    with open(data_file) as f:
        for line in f:
            if '# energy\t' in line:
                de = [ -int(val) for val in line.split()[2:] ]
                break
    return de

de = get_de(data_file)

plt.figure(figsize=fig_size)

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
c = plt.pcolor((de+0.5), e, transitions, vmin=0, vmax=1)
plt.colorbar(c)
lim = min(abs(minde),abs(maxde))
plt.xlim(-lim,lim)
plt.ylim(e.min(), e.max())

plt.xlabel('$\Delta E/\epsilon$')
plt.ylabel('$E/\epsilon$')

plt.tight_layout()
plt.savefig("figs/transitions-example.pdf")
