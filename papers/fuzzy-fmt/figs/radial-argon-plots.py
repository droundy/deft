#!/usr/bin/python3

from __future__ import print_function, division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from pylab import *
import os
import sys
import styles

from cycler import cycler

figscale=1
plt.figure(figsize=(18*figscale, 4.5*figscale))

sigma = 0.3405 #nm     #FIX?  Mikolaj Feb 1967 paper says sigma=3.45 Angstroms = 0.345 nm
sigma_over_R=2**(5/6)

rmaxplot = 4 # upper limit of our plots
#figure(figsize=(18, 4.5))
subplot(1, 3, 1)
title('$n^*=0.8389$, $T^*=0.71$', fontsize=20)
data = loadtxt('figs/YarnellArgon85K.dat')
n =0.02125 # Angstrom^-3
nsig_3 = n*(sigma*10)**3
plot(data[:, 0]/(sigma*10), data[:, 1], label='experiment')

data_dft1 = loadtxt('figs/new-data/radial-lj-0.7100-0.84.dat')
plot(data_dft1[:, 0], data_dft1[:, 1]/nsig_3, label='SFMT')

data_bh1 = loadtxt('figs/new-data/radial-bh-lj-0.7100-0.84.dat')
plot(data_bh1[:, 0], data_bh1[:, 1]/nsig_3, label='BH', ls = ':', color='k')

data_mc1=loadtxt('figs/mc_testp_wca-0.8389-0.7100.dat.gradial')
plot((1/sigma_over_R)*data_mc1[:, 0], data_mc1[:, 1], '--', label='MC')

xlabel(r'$r/\sigma$', fontsize=18)
xlim(0, rmaxplot)
plt.xticks(fontsize=16)
ylabel(r'$g(r)$', fontsize=18)
ylim(0, 3.5)
plt.yticks(fontsize=16)
legend(loc='best', fontsize=18)



subplot(1, 3, 2)
# title('$n^*=0.9570$, $T*=2.48$')
# data2 = loadtxt('figs/EggertArgon0.6GPaX.dat')
# n = 24.23 #atoms/nm^3
# nsig_3 = n*sigma**3
# plot(data2[:,0]/sigma,data2[:,1], label='experiment')
# data_mc2 = loadtxt('figs/mc_testp_wca-0.9570-2.4800.dat.gradial')
# plot((1/sigma_over_R)*data_mc2[:,0], data_mc2[:,1], '--', label='MC')
# data_dft2 = loadtxt('figs/new-data/radial-lj-2.4800-0.96.dat')
# plot(data_dft2[:,0],data_dft2[:,1]/nsig_3, label='DFT')
# xlabel(r'$r/\sigma$')
# xlim(0, rmaxplot)
# legend(loc='best')

title('$n^*=1.0950$, $T^*=2.48$', fontsize=20)
data3 = loadtxt('figs/EggertArgon1.1GPaRAW.dat')
n = 27.74 #atoms/nm^3
nsig_3 = n*sigma**3
plot(data3[:, 0]/sigma, (data3[:, 1]), label='experiment')

data_dft3 = loadtxt('figs/new-data/radial-lj-2.4800-1.09.dat')
plot(data_dft3[:, 0], data_dft3[:, 1]/nsig_3, label='SFMT')

data_bh3 = loadtxt('figs/new-data/radial-bh-lj-2.4800-1.09.dat')
plot(data_bh3[:, 0], data_bh3[:, 1]/nsig_3, label='BH', ls = ':', color='k')

data_mc3 = loadtxt('figs/mc_testp_wca-1.0950-2.4800.dat.gradial')
plot((1/sigma_over_R)*data_mc3[:, 0], data_mc3[:, 1], '--', label='MC')

xlabel(r'$r/\sigma$', fontsize=18)
xlim(0, rmaxplot)
plt.xticks(fontsize=16)
ylabel(r'$g(r)$', fontsize=18)
ylim(0, 3.5)
plt.yticks(fontsize=16)
legend(loc='best', fontsize=18)



subplot(1, 3, 3)
title('$n^*=0.5488$, $T^*=1.235$', fontsize=20)
data4 = loadtxt('figs/Mikolaj-X.dat')
#n = 27.74 #atoms/nm^3
#nsig_3 = n*sigma**3
plot(data4[:, 0]/(sigma*10), data4[:, 1], label='experiment')

data_dft4 = loadtxt('figs/new-data/radial-lj-1.2350-0.58.dat')
plot(data_dft4[:, 0], data_dft4[:, 1]/0.5844, label='SFMT')

data_bh4 = loadtxt('figs/new-data/radial-bh-lj-1.2350-0.58.dat')
plot(data_bh4[:, 0], data_bh4[:, 1]/0.5844, label='BH', ls = ':', color='k')

data_mc5 = loadtxt('figs/mcfcc-0.5844-1.2350.dat.gradial')
plot((1/sigma_over_R)*data_mc5[:, 0], data_mc5[:, 1], '--', label='MC')

# data_dft6 = loadtxt('figs/radial-lj-1.2350-0.58.dat')
# plot((1/sigma_over_R)*data_dft6[:,0],data_dft6[:,1]/.58, '--')

xlabel(r'$r/\sigma$', fontsize=18)
xlim(0, rmaxplot)
plt.xticks(fontsize=16)
ylabel(r'$g(r)$', fontsize=18)
ylim(0, 3.5)
plt.yticks(fontsize=16)
# axvline(1)
# axvline(2**(1.0/6.0))
legend(loc='best', fontsize=18)


tight_layout()
savefig('figs/argon-plots.pdf')

# show()
