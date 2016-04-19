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

all_colors = ['r','g','b','k','c','m','y']

temperature_color = {}
u_temperature_color = {}
eta_color = {}
linetype = {
    1: '-',
    2: '--',
}

Ts = np.arange(0.5, 20.0, 0.25)
Tvals_to_plot = [0.5, 0.75, 1, 1.5, 2, 3, 5, 10]
Tnum_to_plot = [list(Ts).index(T) for T in Tvals_to_plot]

L = 2.84
ww = 1.3
R = 1
i_values = [1,2]
next_color = 0
for i in i_values:
    print 'working on iteration', i
    dbase = 'data/scrunched-ww%04.2f-L%04.2f/i%01d' % (ww, L,  i)
    lndos, energy, Ns = gatherandcalculate.lndos_energy_Ns(dbase)

    V = (L*2**i)**3
    etas = Ns*4*np.pi/3*R**3/V

    Uexc = gatherandcalculate.Uexc(lndos, energy, Ts)
    Fexc = gatherandcalculate.Fexc(dbase, lndos, energy, V, Ts)
    Sexc = gatherandcalculate.Sexc(Uexc, Fexc, Ts)

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

    plt.figure('Uexc-eta')
    for k in Tnum_to_plot:
        ok = Uexc[k,:] != 0
        if Ts[k] in u_temperature_color:
            plt.plot(etas[ok], Uexc[k,:][ok]/Ns[ok], 'x'+linetype[i]+u_temperature_color[Ts[k]])
        else:
            u_temperature_color[Ts[k]] = all_colors[next_color]
            next_color = (next_color+1) % len(all_colors)
            plt.plot(etas[ok], Uexc[k,:][ok]/Ns[ok], linetype[i]+u_temperature_color[Ts[k]],
                     label=r'$T = %g$' % Ts[k])

    plt.figure('HS')
    Sexcs, NNs = gatherandcalculate.Sexc_hardsphere_Ns(dbase)
    print NNs
    plt.plot(NNs*4*np.pi/3*R**3/V, Sexcs/NNs, linetype[i]+'k')


plt.figure('HS')
eta = np.arange(0.01, 0.4, 0.01)
n = eta/(4*np.pi/3*R**3)
Scs = -(4*eta-3*eta**2)/(1-eta)**2
plt.plot(eta, Scs, ':', label='Carnahan-Starling')
plt.title(r'Excess hard-sphere entropy for $L=%g$' % (L))
plt.xlabel(r'$\eta$')
plt.ylabel(r'$S_{exc}/N$')
plt.legend(loc='best')
plt.savefig("figs/Shs-vs-eta.pdf")

plt.figure('Fexc-eta')
plt.title(r'Absolute free energies for $\lambda=%g$, $L=%g$' % (ww,L))
plt.xlim(0,0.6)
plt.ylim(-5,50)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$F_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Fexc-vs-eta.pdf")

plt.figure('Uexc-eta')
plt.title(r'Excess internal energies for $\lambda=%g$, $L=%g$' % (ww,L))
plt.xlabel(r'$\eta$')
plt.ylabel(r'$U_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Uexc-vs-eta.pdf")

plt.figure('Fexc-T')
plt.title(r'Excess free energies for $\lambda=%g$, $L=%g$' % (ww, L))
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$F_{exc}/\epsilon N$')
plt.legend(loc='best')
plt.savefig("figs/Fexc-vs-T.pdf")

