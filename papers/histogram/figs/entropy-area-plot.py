from __future__ import division
import numpy, sys, os, matplotlib
import readnew

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt

if os.path.exists('../data'):
    os.chdir('..')

moviedir = sys.argv[1]
frames = eval(sys.argv[2])

Ss = []
labels = []

for f in frames:
    print f
    e, lndos = readnew.e_lndos('%s/%06d-lnw.dat' % (moviedir, f))
    its = readnew.iterations('%s/%06d-lnw.dat' % (moviedir, f))
    lndos = -lndos
    Smin = lndos[lndos>lndos[0]].min()
    lndos -= Smin
    if len(Ss) > 0:
        Ss.append(lndos - lastdos)
    else:
        Ss = [lndos]
    lastdos = lndos*1
    labels.append('%d iterations' % its)

#    plt.stackplot(e[lndos > lndos[0]], lndos[lndos>lndos[0]])#, label='%d iterations' % its)
plt.stackplot(e, Ss, labels=labels)
plt.plot(e, lndos)

plt.ylim(ymin=0)
plt.legend(loc='best')

plt.show()
