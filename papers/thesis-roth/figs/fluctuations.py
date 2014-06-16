#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')

from pylab import *

from numpy import random
random.seed(0)

ks = arange(0.0, 1.0, .005)
ks[ks>0.5] -= 1

kx, ky = meshgrid(ks, ks)

gaus = random.normal(size = kx.shape)/(kx**2 + ky**2)**.7

for i in range(gaus.shape[0]):
    for j in range(gaus.shape[1]):
        gaus[i,j] = gaus[gaus.shape[0]-i-1,j]
for i in range(gaus.shape[0]):
    for j in range(gaus.shape[1]):
        gaus[i,j] = gaus[gaus.shape[0]-i-1,gaus.shape[1]-j-1]
for i in range(gaus.shape[0]):
    for j in range(gaus.shape[1]):
        gaus[i,j] = gaus[i,gaus.shape[1]-j-1]
        # if i == gaus.shape[0]/2:
        #     gaus[i,j] = 0
        #     print 'i:', i
        # if j == gaus.shape[1]/2:
        #     gaus[i,j] = 0
        #     print 'j:', j

data = ifft2(gaus)

#title('Real space')
print data.shape, gaus.shape
xs = arange(0.0, 1.0, 0.005)
X, Y = meshgrid(xs, xs)
contourf(X, Y, real(data), 60)
xlim(0, 0.5)
ylim(0, 0.5)

# figure()
# title('Fourier space')
# contourf(kx, ky, gaus, 100)

# figure()
# title('kx')
# contourf(kx**2 + ky**2)

# figure()
# title('ky')
# contourf(kx, ky, ky, 100)

xticks([])
yticks([])

tight_layout()

#show()
savefig('figs/fluctuations.pdf')

