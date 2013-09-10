#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if not 'show' in sys.argv:
  matplotlib.use('Agg')
from pylab import *

r = arange(0, 3.0, 0.001)
deltaf = exp(-(r-1.02)**2/.01**2)
deltanear = exp(-(r-1.2)**2/.01**2)
flatf = zeros_like(r)
flatf[r < 2.5] = 1
flatf[r < 1] = 0

sz = (3.0, 2.5)

figure(figsize=sz)
plot(r, deltaf, 'r-')
yticks([])
xticks([0,1,2,3])
xlabel('$r/\sigma$')
ylabel('$\Phi(r)$')
title('contact potential')
tight_layout()
savefig('figs/phi-delta.pdf', transparent=True)

figure(figsize=sz)
plot(r, deltanear, 'r-')
yticks([])
xticks([0,1,2,3])
xlabel('$r/\sigma$')
ylabel('$\Phi(r)$')
title('near contact potential')
tight_layout()
savefig('figs/phi-delta-off.pdf', transparent=True)

figure(figsize=sz)
plot(r, flatf, 'r-')
yticks([])
xticks([0,1,2,3])
ylim(0,1.2)
xlabel('$r/\sigma$')
ylabel('$\Phi(r)$')
title('wide flat potential')
tight_layout()
savefig('figs/phi-flat.pdf', transparent=True)


show()
