#!/usr/bin/python2
from __future__ import division
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import numpy as np

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import gatherandcalculate

def find_common_tangent(x,y):
    los = []
    his = []
    slopes = []
    for i in range(len(x)):
        for j in range(i+1,len(y)):
            slope = (y[j]-y[i])/(x[j]-x[i])
            ynew = y - y[i] - slope*(x-x[i])
            if ynew.min() > -abs(y[i])*1e-5:
                los.append(i)
                his.append(j)
                slopes.append(slope)
    longest = 0
    for i in range(len(los)):
        if his[i]-los[i] > longest:
            lo = los[i]
            hi = his[i]
            longest = x[hi]-x[lo]
            slope = slopes[i]
    return slope

all_colors = ['g','b','r','k','c','m','y']

temperature_color = {}
u_temperature_color = {}
fabs_temperature_color = {}
eta_color = {}
linetype = {
    1: '-',
    2: '--',
    3: 'o:',
}

Ts = np.arange(0.5, 20.0, 0.25)
Tvals_to_plot = [0.5, 0.75, 1, 1.5, 2, 3, 5, 10]
Tnum_to_plot = [list(Ts).index(T) for T in Tvals_to_plot]

L = 2.84
ww = 1.3
R = 1
mu = -20

i_values = [1,2,3]
next_color = 0
next_u_color = 0

for i in i_values:
    print 'working on iteration', i
    dbase = 'data/scrunched-ww%04.2f-L%04.2f/i%01d' % (ww, L,  i)
    lndos, energy, Ns = gatherandcalculate.lndos_energy_Ns(dbase)

    V = (L*2**i)**3
    etas = Ns*4*np.pi/3*R**3/V

    Uexc = gatherandcalculate.Uexc(lndos, energy, Ts)
    Fexc = gatherandcalculate.Fexc(dbase, lndos, energy, V, Ts)
    Sexc = gatherandcalculate.Sexc(Uexc, Fexc, Ts)
    Fabs = gatherandcalculate.Fabs(Fexc, Ts, Ns, V)
    Phiabs = gatherandcalculate.Phiabs(Fabs, mu, Ns)

    plt.figure('Fexc-T')
    for j in range(len(Ns)):
        modulo = 2*8**(i-1)
        if Ns[j] % modulo == 0 and (Fexc[:,j] != 0).any():
            val = Ns[j] // modulo
            c = all_colors[(Ns[j] // modulo) % len(all_colors)]
            if i == i_values[0]:
                plt.plot(Ts, Fexc[:,j]/Ns[j], linetype[i]+c, label=r'$\eta = %.2g$' % etas[j])
            else:
                plt.plot(Ts, Fexc[:,j]/Ns[j], linetype[i]+c)

    plt.figure('Fexc-eta')
    for k in Tnum_to_plot:
        ok = Fexc[k,:] != 0
        if Ts[k] in temperature_color:
            plt.plot(etas[ok], Fexc[k,:][ok]/Ns[ok], linetype[i]+temperature_color[Ts[k]])
        else:
            temperature_color[Ts[k]] = all_colors[next_color]
            next_color = (next_color+1) % len(all_colors)
            plt.plot(etas[ok], Fexc[k,:][ok]/Ns[ok], linetype[i]+temperature_color[Ts[k]],
                     label=r'$T = %g$' % Ts[k])

    plt.figure('Fabs-eta')
    for k in Tnum_to_plot:
        ok = Fabs[k,:] != 0
        if Ts[k] in fabs_temperature_color:
            plt.plot(etas[ok], Fabs[k,:][ok]/V, linetype[i]+fabs_temperature_color[Ts[k]])
        else:
            fabs_temperature_color[Ts[k]] = all_colors[next_color]
            next_color = (next_color+1) % len(all_colors)
            plt.plot(etas[ok], Fabs[k,:][ok]/V, linetype[i]+fabs_temperature_color[Ts[k]],
                     label=r'$T = %g$' % Ts[k])

    plt.figure('Uexc-eta')
    for k in Tnum_to_plot:
        ok = Uexc[k,:] != 0
        if Ts[k] in u_temperature_color:
            plt.plot(etas[ok], Uexc[k,:][ok]/Ns[ok], linetype[i]+u_temperature_color[Ts[k]])
        else:
            u_temperature_color[Ts[k]] = all_colors[next_u_color]
            next_u_color = (next_u_color+1) % len(all_colors)
            plt.plot(etas[ok], Uexc[k,:][ok]/Ns[ok], linetype[i]+u_temperature_color[Ts[k]],
                     label=r'$T = %g$' % Ts[k])

    plt.figure('Phiabs-eta')
    for k in Tnum_to_plot[:3]:
        T = Ts[k]
        mu_here = (Fabs[k,:]/Ns).min()
        mu_here = find_common_tangent(Ns, Fabs[k,:])
        Phi = Fabs[k,:] - mu_here*Ns
        ok = Phi != 0
        if Ts[k] in u_temperature_color:
            plt.plot(etas[ok], Phi[ok]/V, linetype[i]+u_temperature_color[Ts[k]])
        else:
            u_temperature_color[Ts[k]] = all_colors[next_u_color]
            next_u_color = (next_u_color+1) % len(all_colors)
            plt.plot(etas[ok], Phi[ok]/V, linetype[i]+u_temperature_color[Ts[k]],
                     label=r'$T = %g$' % Ts[k])

    plt.figure('Phiabs-T')
    for j in range(len(Ns)):
        modulo = 2*8**(i-1)
        if Ns[j] % modulo == 0 and (Fabs[:,j] != 0).any():
            val = Ns[j] // modulo
            c = all_colors[(Ns[j] // modulo) % len(all_colors)]
            if i == i_values[0]:
                plt.plot(Ts, Phiabs[:,j]/Ns[j], linetype[i]+c, label=r'$\eta = %.2g$' % etas[j])
            else:
                plt.plot(Ts, Phiabs[:,j]/Ns[j], linetype[i]+c)


    plt.figure('HS')
    Sexcs, NNs = gatherandcalculate.Sexc_hardsphere_Ns(dbase)
    print NNs
    plt.plot(NNs*4*np.pi/3*R**3/V, Sexcs/NNs, linetype[i]+'k')


plt.figure('HS')
eta = np.arange(0.01, 0.6, 0.01)
n = eta/(4*np.pi/3*R**3)
Scs = -(4*eta-3*eta**2)/(1-eta)**2
plt.plot(eta, Scs, ':', label='Carnahan-Starling')
plt.title(r'Excess hard-sphere entropy for $L=%g$' % (L))
plt.xlabel(r'$\eta$')
plt.ylabel(r'$S_{exc}/N$')
plt.legend(loc='best')
plt.savefig("figs/Shs-vs-eta.pdf")

plt.figure('Fexc-eta')
plt.title(r'Excess free energies for $\lambda=%g$, $L=%g$' % (ww,L))
plt.xlim(0,0.6)
plt.ylim(-5,50)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$F_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Fexc-vs-eta.pdf")

plt.figure('Fabs-eta')
plt.title(r'Absolute free energies for $\lambda=%g$, $L=%g$' % (ww,L))
plt.xlabel(r'$\eta$')
plt.ylabel(r'$F_{abs}/\epsilon V$')
plt.ylim(-5, 0)
plt.xlim(0,0.6)
plt.legend(loc='best')
plt.savefig("figs/Fabs-vs-eta.pdf")

plt.figure('Uexc-eta')
plt.title(r'Excess internal energies for $\lambda=%g$, $L=%g$' % (ww,L))
plt.xlabel(r'$\eta$')
plt.ylabel(r'$U_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Uexc-vs-eta.pdf")

plt.figure('Fexc-T')
plt.title(r'Excess free energies for $\lambda=%g$, $L=%g$' % (ww, L))
plt.xlim(0,10)
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$F_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Fexc-vs-T.pdf")


plt.figure('Phiabs-eta')
plt.title(r'Absolute grand free energies for $\lambda=%g$, $L=%g$' % (ww, L))
# plt.xlim(0,1)
plt.ylim(-0.02, 0.4)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\Phi_{abs}/\epsilon V$')
plt.legend(loc='best')
plt.savefig("figs/Phiabs-vs-eta.pdf")

# plt.figure('Phiabs-T')
# plt.title(r'Absolute grand free energies for $\lambda=%g$, $L=%g$' % (ww, L))
# plt.xlim(0,10)
# plt.xlabel(r'$kT/\epsilon$')
# plt.ylabel(r'$\Phi_{abs}/\epsilon N$')
# plt.legend(loc='best')
# plt.savefig("figs/Phiabs-vs-T.pdf")
