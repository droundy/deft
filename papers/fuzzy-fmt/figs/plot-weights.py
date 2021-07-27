#!/usr/bin/python3

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import findxi
import scipy.special

T = 1.0
alpha = findxi.find_alpha(T)
Xi = findxi.find_Xi(T)

data = np.loadtxt('figs/weight-functions-%g.dat' % T)

# x	y	z	n	n3	n2	n1	n0	n2x	n2y	n2z	n2xx	n2yy	n2zz	n2xy	n2yz	n2zx
x = data[:,0]
y = data[:,1]
z = data[:,2]
n = data[:,3]
n3 = data[:,4]
n2 = data[:,5]
n1 = data[:,6]
n0 = data[:,7]
n2x = data[:,8]
n2y = data[:,9]
n2z = data[:,10]
n2xx = data[:,11]
n2yy = data[:,12]
n2zz = data[:,13]
n2xy = data[:,14]
n2yz = data[:,15]
n2zx = data[:,16]
print(z)
r = np.sqrt(x**2+y**2+z**2)

def plotone(f,name):
    plt.ylabel(name)
    plt.xlabel('z')
    plt.plot(z, f, label=name)

w2 = np.sqrt(2/np.pi)/Xi*np.exp(-((r-alpha/2)/(Xi/np.sqrt(2)))**2)

plt.figure()
plotone(n, 'n')
plotone(n0, 'n0')
plt.plot(z, w2/(4*np.pi*r**2), '--', label='theory')
plt.plot(z, w2/(4*np.pi*r**2)/2, '--', label='half theory')
plt.ylim(0, 1.2*n0.max())
plt.legend(loc='best')

plt.figure()
plotone(n3, 'n3')
plt.plot(z, 0.5*(1-scipy.special.erf((r-alpha/2)/(Xi/np.sqrt(2)))), '--', label='theory')
plt.legend(loc='best')

plt.figure()
plotone(n2, 'n2')
plt.plot(z, w2, '--', label='theory')
plt.legend(loc='best')

plt.figure()
plotone(n1, 'n1')
plt.plot(z, w2/(4*np.pi*r), '--', label='theory')
plt.legend(loc='best')

plt.figure()
plotone(n2x, 'n2x')
plt.plot(z, w2*x/r, '--', label='x theory')
plotone(n2y, 'n2y')
plt.plot(z, w2*y/r, '--', label='y theory')
plotone(n2z, 'n2z')
plt.plot(z, w2*z/r, '--', label='z theory')
plt.plot(z, w2*z/r/2, '--', label='half z theory')
plt.legend(loc='best')

plt.figure()
plotone(n2xx, 'n2xx')
plt.plot(z, w2*(x**2/r**2 - 1./3), '--', label='xx theory')
plotone(n2yy, 'n2yy')
plt.plot(z, w2*(y**2/r**2 - 1./3), '--', label='yy theory')
plotone(n2zz, 'n2zz')
plt.plot(z, w2*(z**2/r**2 - 1./3), '--', label='zz theory')
plotone(n2xy, 'n2xy')
plotone(n2yz, 'n2yz')
plotone(n2zx, 'n2zx')
plt.legend(loc='best')

plt.show()
