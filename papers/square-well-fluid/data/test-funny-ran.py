from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

dP = 0.000123
P = np.arange(0, 1, dP)

mean = 0.001

b = np.exp(0.5*mean) # np.sqrt(1 + (1.0/mean)**2)
b = 10
invb = 1.0/b

r = (1/(b - P*(b - invb)) - invb)/(b - invb)

print 'mean is', mean
print 'b is', b
print '1/2 invb is', 1/b/2
print 'actual rms is', np.sqrt(sum(r*r)*dP)
print 'actual mean is', sum(r)*dP

print 'mean from formula is', (b*np.log(1/b**2)/(1-b**2) - 1/b)/(b - 1/b)

print 'mean from formula is', (2*b**2*np.log(b)/(b**2-1) - 1)/(b**2 - 1)

plt.plot(P, r)
plt.show()
exit(0)

plt.figure()
bs = np.arange(10.0, 1000.0, 0.1)
plt.loglog(bs, (2*bs**2*np.log(bs)/(bs**2-1) - 1)/(bs**2 - 1), 'b-')

rvals = (2*np.log(bs) - 1)/bs**2
plt.loglog(bs, (2*np.log(bs) - 1)/bs**2, 'r-')
plt.loglog(np.sqrt(np.e + np.e**2*rvals), rvals, 'g--')

plt.show()
