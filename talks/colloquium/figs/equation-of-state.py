#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if not 'show' in sys.argv:
  matplotlib.use('Agg')
from pylab import *

figure(figsize=(7,4))

eos = loadtxt("../../papers/hughes-saft/figs/equation-of-state.dat")
eos_exp = loadtxt("../../papers/hughes-saft/figs/experimental-equation-of-state.dat")

gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

plot(eos[:,2]/gpermL, eos[:,0], 'b-', label='SAFT')
plot(eos[:,3]/gpermL, eos[:,0], 'b-')

plot(eos_exp[:,2]/gpermL, eos_exp[:,0], 'r--', label='experiment')
plot(eos_exp[:,3]/gpermL, eos_exp[:,0], 'r--')

legend(loc='best').draw_frame(False)

ylim(273, 710)
xlim(0, 1.05)
xlabel('$n$ (g/mL)')
ylabel('$T$ (K)')
tight_layout()

savefig('figs/equation-of-state.pdf', transparent=True)
show()
