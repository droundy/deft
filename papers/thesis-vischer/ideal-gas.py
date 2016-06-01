from __future__ import division
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

# Constants

R=3.034E-10 #m, sphere radius for water
M=2.989E-23 #kg, water mass
qe = 1.19E-19

hbar = 1.055E-34 # hbar
kb = 1.381E-23 # boltzmann
N = 60 # spheres
V = N*(4/3)*np.pi*R**3/(.000001) #m^3, volume (size required for \eta = .25) with N=60

# Array generation

colors = ['1', '.75', 'r', '.5', '.25'] # water is red for detail
ms = [9.11E-31, 1.67E-27,M, 80, 1.98E30] # masses to plot
mlabels = ['$m_e$', '$m_P$', '$m_{water}$', '$m_{person}$', '$M_{\odot}$'] # pretty latex labels
Ts = np.arange(1.0,300,.1) # K, temperatures, avoiding 0 because no one will notice
E = np.zeros_like(Ts) 
cv = np.zeros_like(Ts) 

# Calculate free energies for each mass and every temperature

for i in range(0,len(ms)):
    for j in range(0,len(Ts)):
        Lambda = hbar/np.sqrt(kb*ms[i]*Ts[j]/(2*np.pi))
        Z = V/Lambda**3        
        E[j] = (-N*kb*Ts[j]*np.log(Z) + N*kb*Ts[j]*np.log(N) -N*kb*Ts[j])/qe
        cv[j] = E[j]/Ts[j] + 1.5*N*kb
    print "cv is, ", cv
    plt.figure('cv')
    plt.semilogx(Ts, cv/N/kb, color = colors[i], marker = 'o', label = mlabels[i])


plt.figure('cv')        
plt.title('Ideal gas free energy for $N=60$, $\eta = .25$')
plt.ylabel('$E$ (eV)')
plt.xlabel('$T$ (K)')
#plt.ylim(-800, 100)
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig('ideal-F-vs-T.pdf')
plt.show()
