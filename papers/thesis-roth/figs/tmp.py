#!/usr/bin/python

import pylab
import SW
import RG

# Read in data
data = pylab.loadtxt('npart_SW-out.dat')

T = data[:,0]
nvapor = data[:,1]
nliquid = data[:,2]

# Plot the curve
# pylab.subplot(3,1,1)
# pylab.plot(nvapor,T)
# pylab.plot(nliquid,T)
# pylab.text(0.06,3,'Coexistence, SW')
# pylab.ylabel('T')

n = pylab.linspace(0,0.15,10000)
T = 0.07
nparticular = 0.0146971318221

pylab.subplot(2,1,1)
pylab.plot(n,SW.phi(T,n,nparticular),'g-')
pylab.ylabel(r'$\phi_{SW}$')
pylab.text(0.06,.02,'SW; T=%d'%T,color='green')

pylab.subplot(2,1,2)
pylab.plot(n,RG.phi(T,n,nparticular,0),'r-')
pylab.ylabel(r'$\phi_{RG}$')
pylab.text(0.08,.6,'RG; i=0; T=%d'%T,color='red')

pylab.subplot(2,1,2)
pylab.plot(n,RG.phi(T,n,nparticular,1),'b-')
pylab.text(0.08,.4,'RG; i=1; T=%d'%T,color='blue')

# pylab.subplot(2,1,2)
# pylab.plot(n,RG.ftot(T,n,2),'g-')
# pylab.ylabel(r'$f_{RG}$')
# pylab.text(0.08,-2,'RG; i=2; T=%d'%T,color='green')

pylab.xlabel('n')
pylab.show()
