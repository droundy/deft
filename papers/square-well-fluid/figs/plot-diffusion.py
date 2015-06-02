#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readandcompute

ww = float(sys.argv[1])
#arg ww = [1.3]
ffs = eval(sys.argv[2])
#arg ffs = [[0.1,0.2,0.3]]
lenx = float(sys.argv[3])
#arg lenx = [50,100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]

fig, axD = plt.subplots()
axT = plt.twinx()

for ff in ffs:
    basename = 'data/lv/ww%.2f-ff%.2f-%gx%g' % (ww,ff,lenx,lenyz)
    e, diff = readandcompute.e_diffusion_estimate(basename)
    N = readandcompute.read_N(basename);
    axD.plot(e, diff, label=r'$\eta = %g$' % ff)
    axD.axvline(-readandcompute.max_entropy_state(basename)/N, linestyle=':')
    axD.axvline(-readandcompute.min_important_energy(basename)/N, linestyle='--')

    T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
    axT.plot(u/N, T, 'r-')
    axT.set_ylim(0, 3)
    axT.axhline(minT, color='r', linestyle=':')

    e, hist = readandcompute.e_hist(basename)
    axT.plot(e/N, 2.5*hist/hist.max(), 'k:')

axT.spines['right'].set_color('red')
axT.tick_params(axis='y', colors='red')
axT.yaxis.label.set_color('red')

axD.legend()
plt.xlabel(r'$E_{tot}$')
axD.set_ylabel(r'Estimated diffusion constant (actually RMS energy change)')
axD.set_ylim(ymin=0)
axT.set_ylabel(r'$T$')
plt.title(r'RMS $\Delta E$ with $\lambda = %g$' % (ww))

plt.savefig('figs/liquid-vapor-ww%.2f-%gx%g-diffusion.pdf' % (ww,lenx,lenyz))

plt.show()
