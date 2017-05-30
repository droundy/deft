#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readnew

if len(sys.argv) < 5:
    print("Usage: python {} 1.3 0.22 100 10 tmi".format(sys.argv[0]))
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]

if len(sys.argv) > 5:
    method = sys.argv[5]
else:
    method = "-tmi3"
#arg method = ['-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']

plt.figure()


# We manually choose the file with which we wish to compare against.
fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s-movie/000500' % (ww,ff,lenx,lenyz,'-tmi2')

energy_ref,lndos_ref,ps_ref = readnew.e_lndos_ps(fbase)
ps1_ref = numpy.array(ps_ref) + 1

fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s-movie/000500' % (ww,ff,lenx,lenyz,method)

energy,lndos,ps = readnew.e_lndos_ps(fbase)

# Now we subtract to obtain a comparison.

Error = numpy.absolute(lndos - lndos_ref)

plt.semilogx(ps1_ref, Error, label = '$\mid tmi_3 - tmi_2 \mid$')

plt.ylim((0,1.0))
plt.xlim((10,1000000))

plt.xlabel(r'Pessimistic Samples')
plt.ylabel('Error in $\mathcal{D}(\epsilon)$')
plt.title('Error in $\mathcal{D}(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-Error.pdf' % (ww,ff,lenx,lenyz,method))
plt.show()
