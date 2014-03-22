#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 4:
    print 'useage: %s ww ff N' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [200, 1000]

data = numpy.loadtxt("data/periodic-ww%03.1f-ff%04.2f-N%i-E.dat" % (ww, ff, N))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label='DOS ww %.1f ff %.2f N %d' % (ww, ff, N))

        # if args.print: print('  well width:',ww)
        # fig, ax = newFig()
        # for p in [ p for p in paramList
        #            if p.ww == ww and p.ff in args.ff ]:
        #         data = loadtxt(p.Efile,ndmin=2)
        #         
        #         DS = data[:,1]
        #         DS /= sum(DS)
        #         p.max_index = argmax(DS)
        #         semilogy(energy,DS,'.',label=figLabel(p,'N'))
plt.xlabel('$E/N\epsilon$')
plt.ylabel('$D$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-E.pdf" % (ww*10, ff*100, N))

